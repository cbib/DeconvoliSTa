
library('SingleCellExperiment')
library('SpatialExperiment')
library("SpatialDDLS")
library(Seurat)
library(reticulate)

# Default parameter values; overridden by command-line args below.
par <- list(
    annot    = "subclass",
    map_genes = "false"
)
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
par[names(args)] <- args
cat("Given arguments are \n")
print(par)

seurat_obj_scRNA     <- readRDS(par$sc_input)
spatial_seurat_obj   <- readRDS(par$sp_input)
use_gpu   <- par$use_gpu
annot     <- par$annot
map_genes <- par$map_genes

if (use_gpu == "true") {
    use_python("/opt/conda/envs/myenv/bin/python", required = TRUE)
} else {
    install_tf_if_needed <- function() {
        if (!py_module_available("tensorflow")) {
            message("TensorFlow not found. Installing TensorFlow...")
            installTFpython(install.conda = TRUE)
        } else {
            message("TensorFlow is already installed.")
        }
    }
    install_tf_if_needed()
}

py_run_string("
import tensorflow as tf
if tf.test.is_gpu_available():
    print('TensorFlow is using the GPU')
    print('GPU Devices:', tf.config.experimental.list_physical_devices('GPU'))
else:
    print('TensorFlow is not using the GPU')
")

# Check if epochs and batch_size are provided, otherwise set default values
epochs     <- if (is.null(par$epochs))     100L else as.integer(par$epochs)
batch_size <- if (is.null(par$batch_size)) 10L  else as.integer(par$batch_size)

# ---------------------------------------------------------------------------
# Extract raw counts — Seurat v4 (slot=) / v5 (layer=) compatible.
# SeuratObject >= 5.0.0 made slot= defunct; layer= is the v5 API.
# ---------------------------------------------------------------------------
.get_counts <- function(obj) {
    tryCatch(
        GetAssayData(obj, layer = "counts"),
        error = function(e) GetAssayData(obj, slot = "counts")
    )
}

if (inherits(seurat_obj_scRNA, "Seurat")) {
    DefaultAssay(seurat_obj_scRNA) <- "RNA"
    sc_counts <- .get_counts(seurat_obj_scRNA)
} else {
    sc_counts <- seurat_obj_scRNA$counts
}

if (inherits(spatial_seurat_obj, "Seurat")) {
    sp_assay <- names(spatial_seurat_obj@assays)[
        grep("RNA|Spatial", names(spatial_seurat_obj@assays))[1]]
    DefaultAssay(spatial_seurat_obj) <- sp_assay
    sp_counts <- .get_counts(spatial_seurat_obj)
} else {
    sp_counts <- spatial_seurat_obj$counts
}

# ---------------------------------------------------------------------------
# Gene-symbol → Ensembl mapping (map_genes=true, e.g. UKF243 symbols vs GBMap
# Ensembl). Converts the spatial gene axis so it matches the scRNA reference.
# ---------------------------------------------------------------------------
if (map_genes == "true") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
        stop("map_genes=true requires 'org.Hs.eg.db'; use the cbib sp_ddls image.")
    library(org.Hs.eg.db)
    sym2ens   <- mapIds(org.Hs.eg.db,
                        keys    = rownames(sp_counts),
                        keytype = "SYMBOL",
                        column  = "ENSEMBL")
    keep      <- !is.na(sym2ens)
    ens_ids   <- unname(sym2ens[keep])
    sym_names <- names(sym2ens)[keep]
    # Multiple symbols can map to the same Ensembl ID — keep first occurrence.
    uniq_mask <- !duplicated(ens_ids)
    sp_counts <- sp_counts[sym_names[uniq_mask], , drop = FALSE]
    rownames(sp_counts) <- ens_ids[uniq_mask]
    cat("Gene mapping: retained", sum(uniq_mask), "/", length(sym2ens),
        "spatial genes with unique Ensembl IDs\n")
}

# ---------------------------------------------------------------------------
# Subsample scRNA to max_cells_per_type cells per cell type.
# GBMap has 338k cells; converting that to a dense matrix for SpatialDDLS
# requires ~58 GB RAM (OOM). 200 cells/type (~10k total) is standard practice
# and sufficient for training the deconvolution model.
# ---------------------------------------------------------------------------
max_cells_per_type <- 200L
ct_all <- seurat_obj_scRNA@meta.data[[annot]]
set.seed(42L)
keep_idx <- unlist(tapply(seq_along(ct_all), ct_all, function(idx) {
    if (length(idx) > max_cells_per_type) sample(idx, max_cells_per_type) else idx
}), use.names = FALSE)
sc_counts_sub <- sc_counts[, keep_idx, drop = FALSE]
ct_sub        <- ct_all[keep_idx]
cat("scRNA subsampled:", ncol(sc_counts_sub), "cells /",
    length(unique(ct_sub)), "types (max", max_cells_per_type, "per type)\n")

# ---------------------------------------------------------------------------
# Build SpatialExperiment (spatial) and SingleCellExperiment (scRNA) objects.
# Keep sparse matrices (dgCMatrix) — as.matrix() on 338k cells = OOM.
# ---------------------------------------------------------------------------
coordinates <- matrix(0, nrow = ncol(sp_counts), ncol = 2)

sp_input_se_obj <- SpatialExperiment::SpatialExperiment(
    assays        = list(counts = sp_counts),
    rowData       = data.frame(Gene_ID = rownames(sp_counts)),
    colData       = data.frame(Cell_ID = colnames(sp_counts)),
    spatialCoords = coordinates
)

sc_input_sce_obj <- SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = sc_counts_sub),
    colData = data.frame(
        Cell_ID   = colnames(sc_counts_sub),
        Cell_Type = ct_sub
    ),
    rowData = data.frame(Gene_ID = rownames(sc_counts_sub))
)

SDDLS <- createSpatialDDLSobject(
    sc.data           = sc_input_sce_obj,
    sc.cell.ID.column = "Cell_ID",
    sc.gene.ID.column = "Gene_ID",
    st.data           = sp_input_se_obj,
    st.spot.ID.column = "Cell_ID",
    st.gene.ID.column = "Gene_ID",
    project           = "Simul_example",
    sc.filt.genes.cluster = FALSE
)

SDDLS <- genMixedCellProp(
    object           = SDDLS,
    cell.ID.column   = "Cell_ID",
    cell.type.column = "Cell_Type",
    num.sim.spots    = 4000,
    train.freq.cells = 2/3,
    train.freq.spots = 2/3,
    verbose          = TRUE
)

SDDLS <- simMixedProfiles(SDDLS, threads = 8)

SDDLS <- trainDeconvModel(
    object     = SDDLS,
    batch.size = batch_size,
    num.epochs = epochs,
    threads    = 8
)

SDDLS <- loadSTProfiles(
    object            = SDDLS,
    st.data           = sp_input_se_obj,
    st.spot.ID.column = "Cell_ID",
    st.gene.ID.column = "Gene_ID"
)

SDDLS_PROPS <- deconvSpatialDDLS(
    SDDLS, index.st = 1, k.spots = 6, fast.pca = TRUE
)

deconv_matrix <- SDDLS_PROPS@deconv.spots[[1]][6]$Instrinsic.Proportions

# Remove all spaces and dots from cell names, sort them
colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .&-]", "")
deconv_matrix <- deconv_matrix[, sort(colnames(deconv_matrix), method = "shell")]

if (nrow(deconv_matrix) != ncol(sp_counts)) {
    message("The following spots were removed, possibly due to low number of genes: ",
            paste0("'", colnames(sp_counts)[!colnames(sp_counts) %in% rownames(deconv_matrix)], "'",
                   collapse = ", "))
}
rm(sc_counts, sc_counts_sub)  # libère la RAM avant la sérialisation finale

backup_file <- "backup_ddls"
saveRDS(deconv_matrix, file = backup_file)
message("Backup file saved at: ", backup_file)
write.table(deconv_matrix, file = par$output, sep = "\t", quote = FALSE, row.names = FALSE)

 

library('SingleCellExperiment')
library('SpatialExperiment')
library("SpatialDDLS")
library(Seurat)

library(reticulate)

# use_python("/opt/conda/envs/myenv/bin/python", required = TRUE)

# Function to install TensorFlow if not already installed
install_tf_if_needed <- function() {
# Check if TensorFlow is already available
if (!py_module_available("tensorflow")) {
    message("TensorFlow not found. Installing TensorFlow...")
    # Install TensorFlow
    installTFpython(install.conda = TRUE)
} else {
    message("TensorFlow is already installed.")
}
}

# install_tf_if_needed()


# Function to display versions of Python, TensorFlow, and other packages
display_versions <- function() {
python_version <- py_run_string("import sys; print('Python version:', sys.version)")
tensorflow_version <- py_run_string("import tensorflow as tf; print('TensorFlow version:', tf.__version__)")
installed_packages <- py_run_string("import pkg_resources; [print(d) for d in pkg_resources.working_set]")

}

# display_versions()



par <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)
cat("Given arguments are \n")
print(par)
seurat_obj_scRNA <- readRDS(par$sc_input)
spatial_seurat_obj <- readRDS(par$sp_input)

# Check if epochs and batch_size are provided, otherwise set default values
if (is.null(par$epochs)) {
    epochs <- 100 # default value for epochs
} else {
    epochs <- as.integer(par$epochs)
}

if (is.null(par$batch_size)) {
    batch_size <- 10
} else {
    batch_size <- as.integer(par$batch_size)
}


coordinates <- matrix(0, nrow = dim(spatial_seurat_obj$counts)[2], ncol = 2)
counts = spatial_seurat_obj$counts

sp_input_se_obj <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = as.matrix(spatial_seurat_obj$counts)) ,
    rowData = data.frame(Gene_ID = spatial_seurat_obj$counts@Dimnames[[1]]),
    colData = data.frame(Cell_ID = spatial_seurat_obj$counts@Dimnames[[2]]),
    spatialCoords = coordinates
)


sc_input_sce_obj <- SingleCellExperiment::SingleCellExperiment(
    assays = list( counts = as.matrix(seurat_obj_scRNA$RNA@counts)),
    colData = data.frame(
            Cell_ID = seurat_obj_scRNA$RNA@counts@Dimnames[[2]],
            Cell_Type = seurat_obj_scRNA@meta.data$subclass
    ),
    rowData = data.frame(
        Gene_ID = seurat_obj_scRNA$RNA@counts@Dimnames[[1]]
    )
)

SDDLS <- createSpatialDDLSobject(
    sc.data = sc_input_sce_obj,
    sc.cell.ID.column = "Cell_ID",
    sc.gene.ID.column = "Gene_ID",
    st.data = sp_input_se_obj,
    st.spot.ID.column = "Cell_ID",
    st.gene.ID.column = "Gene_ID",
    project = "Simul_example",
    sc.filt.genes.cluster = FALSE
)

#generate mixed cell props profiles from scRNA data for model training
SDDLS <- genMixedCellProp(
    object = SDDLS,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    num.sim.spots = 1000, #5000
    train.freq.cells = 2/3,
    train.freq.spots = 2/3,
    verbose = TRUE
)

SDDLS <- simMixedProfiles(SDDLS, threads = 8)

#Train a deep neural network model using training data from the SpatialDDLS object. This model
#will be used to deconvolute spatial transcriptomics data from the same biological context as the
#single-cell RNA-seq data used to train it
SDDLS <- trainDeconvModel(object = SDDLS,
        batch.size = batch_size,
        num.epochs = epochs,
        threads = 8
    )

#save the model in a HDF5 file  
# saveTrainedModelAsH5(SDDLS, par$model_save_path, overwrite = TRUE)
# model <- loadTrainedModelFromH5(SDDLS, 'ddls_mod_save', reset.slot = FALSE)
# trained.model(SDDLS) <- trained.model(model)

SDDLS <- loadSTProfiles(
    object = SDDLS,
    st.data = sp_input_se_obj,
    st.spot.ID.column = "Cell_ID",
    st.gene.ID.column = "Gene_ID"
)

SDDLS_PROPS <- deconvSpatialDDLS(
    SDDLS, index.st = 1, k.spots = 6, fast.pca = TRUE
)


deconv_matrix <- SDDLS_PROPS@deconv.spots[[1]][6]$Instrinsic.Proportions

# Remove all spaces and dots from cell names, sort them
colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .&-]", "")
deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

if (nrow(deconv_matrix) != ncol(spatial_seurat_obj$counts)){
    message("The following rows were removed, possibly due to low number of genes: ",
            paste0("'", colnames(spatial_seurat_obj$counts)[!colnames(spatial_seurat_obj$counts) %in% rownames(deconv_matrix)], "'", collapse=", "))
}




# backup_file <- tempfile(fileext = ".rds")
backup_file <- "backup_ddls"
saveRDS(deconv_matrix, file = backup_file)
message("Backup file saved at: ", backup_file)
write.table(deconv_matrix, file=par$output, sep="\t", quote=FALSE, row.names=FALSE)




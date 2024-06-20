# if (!requireNamespace("SpatialDDLSdata", quietly = TRUE)) {
#     remotes::install_github('diegommcc/SpatialDDLSdata')
# }
# library("SpatialDDLSdata")
# SingleCellExperiment with scRNA-seq

# data(MouseDLN.SCE) 
# # SpatialExperiment with spatial transcriptomics data
# data(MouseDLN.ST)

# if (!requireNamespace('R.utils', quietly = TRUE)) {
    
# }

library('SingleCellExperiment')
library('SpatialExperiment')
library("SpatialDDLS")

par <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)
print(par)
library(Seurat)
seurat_obj_scRNA <- readRDS(par$sc_input)

sc_input_sce_obj =  as.SingleCellExperiment(seurat_obj_scRNA )



sc_input_ddls <- createSpatialDDLSobject(
  sc.data = par$sc_input, 
  sc.cell.ID.column = "CellID", 
  sc.gene.ID.column = "GeneSymbol",
  sc.cell.type.column = "broad_cell_types",
  st.data = par$sp_input,
  st.spot.ID.column = "CellID",
  st.gene.ID.column = "GeneSymbol",
  sc.filt.genes.cluster = TRUE, 
  sc.n.genes.per.cluster = 150,
  sc.min.mean.counts = 2
)



set.seed(123)
sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(      
        counts = matrix(
            rpois(30, lambda = 5), nrow = 15, ncol = 20,
                dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(20))))
            ),
    colData = data.frame(
            Cell_ID = paste0("RHC", seq(20)),
            Cell_Type = sample(x = paste0("CellType", seq(6)), size = 20,
            replace = TRUE)
    ),
    rowData = data.frame(
        Gene_ID = paste0("Gene", seq(15))
    )
)
SDDLS <- createSpatialDDLSobject(
    sc.data = sce,
    sc.cell.ID.column = "Cell_ID",
    sc.gene.ID.column = "Gene_ID",
    sc.filt.genes.cluster = FALSE
)
SDDLS <- genMixedCellProp(
    object = SDDLS,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    num.sim.spots = 50,
    train.freq.cells = 2/3,
    train.freq.spots = 2/3,
    verbose = TRUE
)
SDDLS <- simMixedProfiles(SDDLS)
    # training of SDDLS model
    SDDLS <- trainDeconvModel(object = SDDLS,
    batch.size = 15,
    num.epochs = 5
)
# simulating spatial data
ngenes <- sample(3:40, size = 1)
ncells <- sample(10:40, size = 1)
counts <- matrix(
rpois(ngenes * ncells, lambda = 5), ncol = ncells,
    dimnames = list(paste0("Gene", seq(ngenes)), paste0("Spot", seq(ncells)))
)
coordinates <- matrix(
    rep(c(1, 2), ncells), ncol = 2
)
st <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = as.matrix(counts)),
    rowData = data.frame(Gene_ID = paste0("Gene", seq(ngenes))),
    colData = data.frame(Cell_ID = paste0("Spot", seq(ncells))),
    spatialCoords = coordinates
)
SDDLS <- loadSTProfiles(
    object = SDDLS,
    st.data = st,
    st.spot.ID.column = "Cell_ID",
    st.gene.ID.column = "Gene_ID"
)
# simplify arguments
simplify <- list(CellGroup1 = c("CellType1", "CellType2", "CellType4"),
    CellGroup2 = c("CellType3", "CellType5"))
SDDLS <- deconvSpatialDDLS(
    object = SDDLS,
    index.st = 1,
    simplify.set = simplify,
    simplify.majority = simplify
)
# mouseDLN.SDDLS <- genMixedCellProp(
#   mouseDLN.SDDLS,
#   cell.ID.column = "CellID",
#   cell.type.column = "broad_cell_types",
#   num.sim.spots = 10000,
#   n.cells = 50,
#   min.zero.prop = 5,
#   balanced.type.cells = TRUE
# )


# mouseDLN.SDDLS <- simMixedProfiles(mouseDLN.SDDLS, threads = 3)
# mouseDLN.SDDLS <- trainDeconvModel(
#   mouseDLN.SDDLS,
#   verbose = FALSE
# ) 






# mouseDLN.SDDLS <- deconvSpatialDDLS(
#   mouseDLN.SDDLS, index.st = 1, k.spots = 6, fast.pca = TRUE
# )
# plotSpatialPropAll(mouseDLN.SDDLS, index.st = 1)


# remotes::install_version('Matrix', version = '1.6.0')
# remotes::install_cran(c('SpatialExperiment', 'SingleCellExperiment', 'SummarizedExperiment', 'zinbwave', 'S4Vectors', 'scran', 'scuttle'))
# remotes::install_github('diegommcc/SpatialDDLS', upgrade='never')
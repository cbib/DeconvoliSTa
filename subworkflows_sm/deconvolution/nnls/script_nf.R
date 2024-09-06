#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(Seurat)
library(nnls)
library(magrittr)

# library(org.Hs.eg.db)
# convert_query_geneSymbol_to_ensemblID <- function(spatial_query){
  
#   ## Function to convert gene symbols in query to ensemble IDs [e.g. GBMap is in ENSEMBL while Visium is in Symbol]
#   ## The following code chunk is to be applied only if the reference single cell dataset
#   ## contains ENSEMBL identifiers instead of gene symbols (as in Visium data)
  
#   master_gene_table <- mapIds(org.Hs.eg.db, keys = spatial_query@Dimnames[[1]], keytype = "SYMBOL", column="ENSEMBL")
#   #master_gene_table <- as.data.frame(master_gene_table)
#     cat ("oaak\n")
#   # get the Ensembl ids with gene symbols i.e. remove those with NA's for gene symbols
#   inds <- which(!is.na(master_gene_table))
#   found_genes <- master_gene_table[inds]
  
#   # subset your data frame based on the found_genes
#   # df2 <- spatial_query@counts[names(found_genes), ]
#   spatial_query@Dimnames[[1]] <- found_genes
  
#   return(spatial_query)
# }


# Accept command line arguments
par <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

print(par)

## START ##
# READING IN INPUT
cat("Reading input scRNA-seq reference from", par$sc_input, "\n")
seurat_obj_scRNA <- readRDS(par$sc_input)
DefaultAssay(seurat_obj_scRNA) <- "RNA"
ncelltypes <- length(unique(seurat_obj_scRNA[[par$annot, drop=TRUE]]))
cat("Found ", ncelltypes, "cell types in the reference.\n")

cat("Reading input spatial data from", par$sp_input, "\n")
spatial_data <- readRDS(par$sp_input)

cat("Getting count matrix of spatial data...\n")
if (class(spatial_data) != "Seurat"){
  spatial_data <- spatial_data$counts
} else {
  DefaultAssay(spatial_data) <- names(spatial_data@assays)[grep("RNA|Spatial",names(spatial_data@assays))[1]]
  spatial_data <- GetAssayData(spatial_data, slot="counts")
}

# spatial_data = convert_query_geneSymbol_to_ensemblID(spatial_data)

# RUNNING THE METHOD #
cat("Calculating basis matrix based on mean average expression per cell type...\n")
basis_matrix <- lapply(SplitObject(seurat_obj_scRNA, split.by=par$annot),
    function(seurat_obj_scRNA_ct) {
        GetAssayData(seurat_obj_scRNA_ct) %>% rowMeans
}) %>% do.call(cbind, .)

cat("Only selecting overlapping genes...\n")
intersect_genes <- intersect(rownames(basis_matrix), rownames(spatial_data))
basis_matrix <- basis_matrix[intersect_genes,]
spatial_data <- spatial_data[intersect_genes,]

cat("Running NNLS...\n")
deconv_matrix <- apply(spatial_data, 2, function(each_spot) {
  x <- nnls(basis_matrix, each_spot)$x
  x/sum(x)
}) %>% t %>% set_colnames(colnames(basis_matrix))

# PRINTING RESULTS #
cat("Printing results...\n")
# Remove all spaces and dots from cell names, sort them
colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .&-]", "")
deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

write.table(deconv_matrix, file=par$output, sep="\t", quote=FALSE, row.names=TRUE)
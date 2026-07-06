#!/usr/bin/env Rscript
# Standard (graph-based, transcriptomic) Seurat clustering of a Visium sample, used alongside the
# spatially-aware BayesSpace clustering in the visualization (the Seurat <-> BayesSpace toggle:
# transcriptomic clusters vs spatial domains).
#
# Output: a CSV with a `BayesSpace` column (cluster per spot, barcode as index) -- the same format
# the visualizer expects for either clustering (the column name is reused for both).
#
# Usage:
#   Rscript seurat_cluster.R --sp_input sample.rds --output seurat_clusters.csv [--resolution 0.8] [--n_pcs 20]

suppressMessages(library(Seurat))
set.seed(42)

par <- list(sp_input = NA, output = "seurat_clusters.csv", resolution = 0.8, n_pcs = 20)
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
par[names(args)] <- args
par$resolution <- as.numeric(par$resolution)
par$n_pcs <- as.integer(par$n_pcs)

cat("Reading spatial object from", par$sp_input, "...\n")
obj <- readRDS(par$sp_input)
DefaultAssay(obj) <- names(obj@assays)[grep("Spatial|RNA", names(obj@assays))[1]]

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = par$n_pcs, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:par$n_pcs, verbose = FALSE)
obj <- FindClusters(obj, resolution = par$resolution, verbose = FALSE)

# Cluster per spot, barcode as index, column named `BayesSpace`. Seurat labels clusters 0-based;
# shift to 1-based so it matches BayesSpace's convention (no cluster 0, which the visualizer paints
# grey #CCCCCC as a neutral color) -- both clusterings then start at 1.
df <- data.frame(BayesSpace = as.integer(as.character(Idents(obj))) + 1L)
rownames(df) <- colnames(obj)
write.csv(df, par$output)
cat("Wrote", par$output, ":", nrow(df), "spots,", length(unique(df$BayesSpace)), "clusters\n")

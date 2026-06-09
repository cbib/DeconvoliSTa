#!/usr/bin/env Rscript
# BayesSpace spatial clustering of a Visium sample (Seurat .rds object)
# -> CSV with a `BayesSpace` column (cluster per spot), the format expected by
#    sp_visualizer.py (data_clustered input).
#
# Usage:
#   Rscript bayesspace_cluster.R --sp_input UKF243_T_ST_1_raw.rds --output seurat_metadata.csv \
#                                [--q auto|<int>] [--qs_min 3] [--qs_max 12] [--nrep 10000]
#
# q="auto" (default): qTune over [qs_min, qs_max] + elbow pick (kneedle). Override with an integer.

suppressMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(BayesSpace)
})
set.seed(42)

par <- list(sp_input = NA, output = "seurat_metadata.csv",
            q = "auto", qs_min = 3, qs_max = 12, nrep = 10000, burn_in = 1000,
            n_pcs = 15, n_hvgs = 2000, cores = 1)
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
par[names(args)] <- args
for (k in c("qs_min","qs_max","nrep","burn_in","n_pcs","n_hvgs","cores")) par[[k]] <- as.integer(par[[k]])

# --- Read the spatial object + build the SingleCellExperiment ---
cat("Reading spatial object from", par$sp_input, "...\n")
obj <- readRDS(par$sp_input)
counts <- GetAssayData(obj, slot = "counts")
co <- obj@images[[1]]@coordinates    # tissue, row, col, imagerow, imagecol
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = DataFrame(row = co$row, col = co$col, row.names = colnames(obj))
)
cat("SCE:", ncol(sce), "spots,", nrow(sce), "genes\n")

cat("spatialPreprocess (PCA, HVGs)...\n")
sce <- spatialPreprocess(sce, platform = "Visium",
                         n.PCs = par$n_pcs, n.HVGs = par$n_hvgs, log.normalize = TRUE)

# --- Choose q ---
pick_elbow <- function(q, ll) {
  # kneedle: point of the (q, loglik) curve farthest from the line joining the extremes
  x <- (q - min(q)) / (max(q) - min(q))
  y <- (ll - min(ll)) / (max(ll) - min(ll))
  x1 <- x[1]; y1 <- y[1]; x2 <- x[length(x)]; y2 <- y[length(y)]
  d <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) / sqrt((y2 - y1)^2 + (x2 - x1)^2)
  q[which.max(d)]
}

if (identical(tolower(as.character(par$q)), "auto")) {
  qs <- seq(par$qs_min, par$qs_max)
  cat("qTune over q =", min(qs), "..", max(qs), "(may be slow)...\n")
  sce <- qTune(sce, qs = qs, platform = "Visium", burn.in = par$burn_in, nrep = par$nrep, cores = par$cores)
  ll <- attr(sce, "q.logliks")
  q_chosen <- pick_elbow(ll$q, ll$loglik)
  cat("q chosen (elbow) =", q_chosen, "\n")
  print(ll)
} else {
  q_chosen <- as.integer(par$q)
  cat("q forced =", q_chosen, "\n")
}

cat("spatialCluster (q =", q_chosen, ", nrep =", par$nrep, ")...\n")
sce <- spatialCluster(sce, q = q_chosen, platform = "Visium",
                      nrep = par$nrep, burn.in = par$burn_in)

# --- Output in the format expected by the visu (BayesSpace column, barcode as index) ---
df <- data.frame(BayesSpace = as.integer(sce$spatial.cluster))
rownames(df) <- colnames(sce)
write.csv(df, par$output)
cat("Wrote", par$output, ":", nrow(df), "spots,", length(unique(df$BayesSpace)), "clusters\n")

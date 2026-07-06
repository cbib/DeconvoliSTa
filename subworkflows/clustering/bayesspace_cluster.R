#!/usr/bin/env Rscript
# BayesSpace spatial clustering of a Visium sample (Seurat .rds object)
# -> CSV with a `BayesSpace` column (cluster per spot), the format expected by
#    sp_visualizer.py (data_clustered input).
#
# Usage:
#   Rscript bayesspace_cluster.R --sp_input UKF243_T_ST_1_raw.rds --output seurat_metadata.csv \
#                                [--q auto|<int>|<int,int,...>] [--qs_min 3] [--qs_max 12] [--nrep 10000]
#
# --q "auto" (default): qTune over [qs_min, qs_max], print a per-step %-improvement table, and
#   pick a kneedle elbow (heuristic only). --q 9: cluster directly at q=9. --q "5,9": cluster at
#   both, writing one file per q (suffixed _q5 / _q9) for side-by-side comparison in the visu.

suppressMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(BayesSpace)
  library(ggplot2)
})
set.seed(42)

par <- list(sp_input = NA, output = "seurat_metadata.csv",
            q = "auto", qs_min = 3, qs_max = 12, nrep = 10000, burn_in = 1000,
            n_pcs = 15, n_hvgs = 2000)
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
par[names(args)] <- args
for (k in c("qs_min","qs_max","nrep","burn_in","n_pcs","n_hvgs")) par[[k]] <- as.integer(par[[k]])

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

# Parse --q: "auto" (qTune + kneedle), or one/several integers (comma-separated) to cluster
# directly. Several integers -> one output file per q (suffixed _q<val>), so different cluster
# counts can be compared side by side in the visu.
q_arg <- tolower(as.character(par$q))
if (identical(q_arg, "auto")) {
  qs <- seq(par$qs_min, par$qs_max)
  cat("qTune over q =", min(qs), "..", max(qs), "(may be slow)...\n")
  # qTune (BayesSpace 1.4.1) has no parallelism: it fits one model per q sequentially.
  # `platform` is forwarded to spatialCluster via `...`.
  sce <- qTune(sce, qs = qs, platform = "Visium", burn.in = par$burn_in, nrep = par$nrep)
  ll <- attr(sce, "q.logliks")

  # Per-step % improvement of the (pseudo) log-likelihood, the way our collaborators inspect it:
  #   pct(q) = |loglik(q) - loglik(q-1)| / |loglik(q-1)| * 100
  # This table is the REAL decision aid: the kneedle elbow below is only a heuristic and tends to
  # UNDERESTIMATE q on the smooth, never-plateauing curves typical of heterogeneous tissue
  # (glioblastoma). Pick the q where the extra cluster stops being biologically worth it.
  pct <- abs(diff(ll$loglik)) / abs(head(ll$loglik, -1)) * 100
  imp_tbl <- data.frame(q_from = head(ll$q, -1), q_to = ll$q[-1],
                        loglik = round(ll$loglik[-1], 2), pct_improvement = round(pct, 2))
  cat("Log-likelihood improvement per added cluster:\n")
  print(imp_tbl, row.names = FALSE)

  q_chosen <- pick_elbow(ll$q, ll$loglik)
  cat("kneedle elbow suggests q =", q_chosen,
      "(heuristic only -- cross-check with the table above and the qPlot)\n")

  # Save the qPlot so the chosen q can be sanity-checked by eye, then override with --q <int>.
  qplot_path <- sub("\\.csv$", "_qplot.png", par$output)
  tryCatch({
    ggsave(qplot_path, plot = qPlot(sce) + geom_vline(xintercept = q_chosen, linetype = "dashed"),
           width = 6, height = 4, dpi = 120)
    cat("qPlot saved to", qplot_path, "\n")
  }, error = function(e) cat("Could not save qPlot:", conditionMessage(e), "\n"))

  qs_to_run <- q_chosen
} else {
  qs_to_run <- as.integer(strsplit(q_arg, ",")[[1]])
  cat("q forced =", paste(qs_to_run, collapse = ", "), "\n")
}

# --- Spatial clustering: one run (and one output file) per requested q ---
# Output format expected by the visu: a `BayesSpace` column with the spot barcode as index.
multi <- length(qs_to_run) > 1
for (qv in qs_to_run) {
  cat("spatialCluster (q =", qv, ", nrep =", par$nrep, ")...\n")
  sce_c <- spatialCluster(sce, q = qv, platform = "Visium",
                          nrep = par$nrep, burn.in = par$burn_in)
  df <- data.frame(BayesSpace = as.integer(sce_c$spatial.cluster))
  rownames(df) <- colnames(sce_c)
  out <- if (multi) sub("\\.csv$", paste0("_q", qv, ".csv"), par$output) else par$output
  write.csv(df, out)
  cat("Wrote", out, ":", nrow(df), "spots,", length(unique(df$BayesSpace)), "clusters\n")
}

# BayesSpace spatial-clustering image (bayesspace_cluster.R + seurat_cluster.R +
# extract_visium_assets.R). Docker equivalent of bayesspace.def -> published to
# ghcr.io/cbib/deconvolista-bayesspace by .github/workflows/build-images.yml.
#
# Overlay on the Seurat image (R 4.1 / Bioconductor 3.14) to add BayesSpace. The R
# scripts are mounted at runtime (Singularity --bind), so the image only ships deps.
FROM csangara/sp_spatialdwls:latest

# Auto-link the published package to the repo on GHCR (provenance + visibility).
LABEL org.opencontainers.image.source=https://github.com/cbib/DeconvoliSTa

# R 4.1.3 in the base image: the latest CRAN versions require R >= 4.3 (e.g. xgboost).
# -> point CRAN at a dated Posit snapshot (2022-10-01) = R 4.1-compatible versions.
# BiocManager installs BayesSpace (Bioc 3.14) + its CRAN deps (xgboost...) from the snapshot.
RUN R -e "options(repos=c(CRAN='https://packagemanager.posit.co/cran/2022-10-01')); \
    if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); \
    BiocManager::install('BayesSpace', version='3.14', ask=FALSE, update=FALSE)"

# Sanity check at build time (fails the build if BayesSpace didn't install).
RUN R -e "library(BayesSpace); library(Seurat); cat('BayesSpace', as.character(packageVersion('BayesSpace')), '+ Seurat OK\n')"

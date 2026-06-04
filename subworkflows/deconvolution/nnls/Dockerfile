FROM csangara/seurat:4.1.0
RUN R -e "install.packages(c('nnls'))"
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); BiocManager::install('org.Hs.eg.db', ask=FALSE)"
# FROM satijalab/seurat as seurat_image

# # Créer un répertoire pour les bibliothèques R
# RUN mkdir /R_libs

# # Copier les bibliothèques R de l'image Seurat
# RUN cp -r /usr/local/lib/R/site-library/* /R_libs/
# FROM r-base:4.4.0
# # Copier les bibliothèques R de l'étape précédente
# COPY --from=seurat_image /R_libs /usr/local/lib/R/site-library/

# RUN ls -l /usr/local/lib/R/site-library/ 
# RUN if [ ! -d "/usr/local/lib/R/site-library/Seurat" ]; then \
#         echo "Seurat package not found. Exiting build." && exit 1; \
#     fi
# RUN apt-get update && apt-get install libxml2 libglpk-dev libicu-dev  -y 
# RUN wget http://archive.ubuntu.com/ubuntu/pool/main/i/icu/libicu66_66.1-2ubuntu2_amd64.deb
# RUN  dpkg -i libicu66_66.1-2ubuntu2_amd64.deb

# # RUN R -e "install.packages('Seurat')"
# RUN R -e "library('Seurat')"


# # Install R dependencies
# # RUN install2.r --error BiocManager version = '3.19'
# RUN R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
# run r -e "BiocManager::install(version = '3.19', ask = FALSE)"
# RUN R -e "library('BiocManager')"

# # RUN R -e "BiocManager::install(c('SpatialExperiment', 'SingleCellExperiment', 'SummarizedExperiment', 'zinbwave', 'S4Vectors', 'scran', 'scuttle'))"
# RUN R -e "BiocManager::install('BiocGenerics')"
# RUN R -e "library('BiocGenerics')"

# RUN R -e  "BiocManager::install('S4Vectors')"
# RUN R -e "library('S4Vectors')"

# RUN R -e "install.packages('crayon')"
# RUN R -e "library('crayon')"

# RUN R -e "install.packages('abind')"
# RUN R -e "library('abind')"

# RUN R -e  "BiocManager::install('S4Arrays')"
# RUN R -e "library('S4Arrays')"

# RUN R -e  "BiocManager::install('Biobase')"
# RUN R -e "library('Biobase')"

# RUN R -e  "BiocManager::install('GenomeInfoDbData')"
# RUN R -e "library('GenomeInfoDbData')"
# RUN apt-get update && apt-get install libcurl4-openssl-dev -y
# RUN  apt-get update && apt-get install libssl-dev -y
# RUN R -e "install.packages('httr')"
# RUN R -e "library('httr')"

# RUN R -e "install.packages('jsonlite')"
# RUN R -e "library('jsonlite')"


# RUN R -e  "BiocManager::install('UCSC.utils')"
# RUN R -e "library('UCSC.utils')"


# RUN R -e  "BiocManager::install('GenomeInfoDb')"
# RUN R -e "library('GenomeInfoDb')"

# RUN R -e  "BiocManager::install('GenomicRanges')"
# RUN R -e "library('GenomicRanges')"

# RUN R -e  "BiocManager::install('MatrixGenerics')"
# RUN R -e "library('MatrixGenerics')"

# RUN R -e  "BiocManager::install('SummarizedExperiment')"
# RUN R -e "library('SummarizedExperiment')"

# RUN R -e  "BiocManager::install('BiocFileCache')"
# RUN R -e "library('BiocFileCache')"


# RUN R -e  "BiocManager::install('SingleCellExperiment')"
# RUN R -e "library('SingleCellExperiment')"


# # RUN apt-get update && apt-get install libglib2.0-0t64=2.80.2-2  libglib2.0-bin=2.80.2-2 libmagick++-dev -y --allow-downgrades
# RUN apt-get update && apt-get install libglib2.0-0t64=2.80.2-2  libglib2.0-bin=2.80.2-2  -y --allow-downgrades

# RUN apt-get update && apt-get install libmagick++-dev  -y --allow-downgrades

# # RUN R -e "if ('SpatialExperiment' %in% rownames(installed.packages())) remove.packages('SpatialExperiment')"

# RUN R -e "BiocManager::install('SpatialExperiment', ask=FALSE)"
# RUN R -e "library('SpatialExperiment')"


# RUN R -e "install.packages('remotes')"
# RUN R -e "remotes::install_version('Matrix', version = '1.6.0')"

# RUN R -e "BiocManager::install('AnnotationDbi')"
# RUN R -e "library('AnnotationDbi')"

# RUN R -e "BiocManager::install('edgeR')"
# RUN R -e "library('edgeR')"

# RUN R -e "BiocManager::install('genefilter')"
# RUN R -e "library('genefilter')"

# RUN R -e "BiocManager::install('BiocParallel')"
# RUN R -e "library('BiocParallel')"

# RUN R -e "remotes::install_version('Matrix')"

# RUN R -e "BiocManager::install('zinbwave')"
# RUN R -e "library('zinbwave')"

# #s  cran dep
# RUN R -e "BiocManager::install('metapod')"
# RUN R -e "library('metapod')"

# RUN R -e "BiocManager::install('BiocNeighbors')"
# RUN R -e "library('BiocNeighbors')"

# RUN R -e "BiocManager::install('bluster')"
# RUN R -e "library('bluster')"

# RUN R -e "BiocManager::install('beachmat')"
# RUN R -e "library('beachmat')"

# RUN R -e "BiocManager::install('ScaledMatrix')"
# RUN R -e "library('ScaledMatrix')"

# RUN R -e "BiocManager::install('DelayedMatrixStats')"
# RUN R -e "library('DelayedMatrixStats')"

# RUN R -e "BiocManager::install('BiocSingular')"
# RUN R -e "library('BiocSingular')"

# RUN R -e "BiocManager::install('limma')"
# RUN R -e "library('limma')"

# RUN R -e "BiocManager::install('scran')"
# RUN R -e "library('scran')"

# RUN apt-get update && apt-get install cmake -y
# RUN R -e "BiocManager::install('SpatialExperiment')"
# RUN R -e "library('SpatialExperiment')"

# RUN R -e "install.packages('SpatialDDLS')"
# RUN R -e "library('SpatialDDLS')"

# RUN R -e  "install.packages('R.utils')"
# RUN R -e "library('R.utils')"

# # Download and install Python 3.8.0
# RUN wget https://www.python.org/ftp/python/3.8.0/Python-3.8.0.tgz \
#     && tar -xf Python-3.8.0.tgz \
#     && cd Python-3.8.0 \
#     && ./configure --enable-optimizations \
#     && make -j 8 \
#     && make altinstall \
#     && cd .. \
#     && rm -rf Python-3.8.0 Python-3.8.0.tgz

# # Update alternatives to use python3.8 as the default python
# RUN update-alternatives --install /usr/bin/python python /usr/local/bin/python3.8 1

# # Install pip for the newly installed Python version
# RUN wget https://bootstrap.pypa.io/get-pip.py \
#     && /usr/local/bin/python3.8 get-pip.py

# # Create and activate a virtual environment using Python 3.8
# RUN /usr/local/bin/python3.8 -m venv /venv
# ENV PATH="/venv/bin:$PATH"


# RUN apt-get update && apt-get install  pkg-config libhdf5-dev -y


# # Install TensorFlow within the virtual environment
# RUN pip install tensorflow








# Start with NVIDIA CUDA base image
FROM nvidia/cuda:11.8.0-base-ubuntu22.04

# Install necessary system dependencies
RUN apt-get update && apt-get install -y \
    software-properties-common \
    dirmngr \
    wget \
    gpg-agent \
    && rm -rf /var/lib/apt/lists/*

# Add R repository
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN apt-get update && apt-get install libxml2 libglpk-dev libicu-dev tzdata -y 
RUN wget http://archive.ubuntu.com/ubuntu/pool/main/i/icu/libicu66_66.1-2ubuntu2_amd64.deb
RUN  dpkg -i libicu66_66.1-2ubuntu2_amd64.deb


# Add the CRAN repository and key
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

# Install R 4.4.0
RUN apt-get update && apt-get install -y r-base=4.4.0-1.2004.0 r-recommended=4.4.0-1.2004.0

# Copy R libraries from the previous image
COPY --from=satijalab/seurat /usr/local/lib/R/site-library /usr/local/lib/R/site-library/

# Install CUDA-specific packages if needed
RUN apt-get update && apt-get install -y cuda-toolkit-11-8


# RUN R -e "install.packages('Seurat')"
RUN R -e "library('Seurat')"

# Install NVIDIA Container Toolkit
RUN apt-get update && apt-get install -y nvidia-container-toolkit

# Set up NVIDIA runtime
ENV NVIDIA_VISIBLE_DEVICES all
ENV NVIDIA_DRIVER_CAPABILITIES compute,utility



RUN apt-get update && apt-get install libffi-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev -y
# Download and install Python 3.8.0
RUN wget https://www.python.org/ftp/python/3.8.0/Python-3.8.0.tgz \
    && tar -xf Python-3.8.0.tgz \
    && cd Python-3.8.0 \
    && ./configure --enable-optimizations \
    && make -j 8 \
    && make altinstall \
    && cd .. \
    && rm -rf Python-3.8.0 Python-3.8.0.tgz

# Update alternatives to use python3.8 as the default python
RUN update-alternatives --install /usr/bin/python python /usr/local/bin/python3.8 1
RUN apt-get update && apt-get install libcurl4-openssl-dev -y
RUN  apt-get update && apt-get install libssl-dev -y


# Install pip for the newly installed Python version
RUN wget https://bootstrap.pypa.io/get-pip.py \
    && /usr/local/bin/python3.8 get-pip.py

# Create and activate a virtual environment using Python 3.8
RUN /usr/local/bin/python3.8 -m venv /venv
ENV PATH="/venv/bin:$PATH"


RUN apt-get update && apt-get install  pkg-config libhdf5-dev -y


# Install TensorFlow within the virtual environment
RUN pip install tensorflow




# Test if GPU is available
RUN python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"

# Install R dependencies
# RUN install2.r --error BiocManager version = '3.19'
# RUN R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
run R -e "BiocManager::install(version = '3.19', ask = FALSE)"
RUN R -e "library('BiocManager')"

# RUN R -e "BiocManager::install(c('SpatialExperiment', 'SingleCellExperiment', 'SummarizedExperiment', 'zinbwave', 'S4Vectors', 'scran', 'scuttle'))"
RUN R -e "BiocManager::install('BiocGenerics')"
RUN R -e "library('BiocGenerics')"

RUN R -e  "BiocManager::install('S4Vectors')"
RUN R -e "library('S4Vectors')"

RUN R -e "install.packages('crayon')"
RUN R -e "library('crayon')"

RUN R -e "install.packages('abind')"
RUN R -e "library('abind')"

RUN R -e  "BiocManager::install('S4Arrays')"
RUN R -e "library('S4Arrays')"

RUN R -e  "BiocManager::install('Biobase')"
RUN R -e "library('Biobase')"

RUN R -e  "BiocManager::install('GenomeInfoDbData')"
RUN R -e "library('GenomeInfoDbData')"

RUN R -e "install.packages('httr')"
RUN R -e "library('httr')"

RUN R -e "install.packages('jsonlite')"
RUN R -e "library('jsonlite')"


RUN R -e  "BiocManager::install('UCSC.utils')"
RUN R -e "library('UCSC.utils')"


RUN R -e  "BiocManager::install('GenomeInfoDb')"
RUN R -e "library('GenomeInfoDb')"

RUN R -e  "BiocManager::install('GenomicRanges')"
RUN R -e "library('GenomicRanges')"

RUN R -e  "BiocManager::install('MatrixGenerics')"
RUN R -e "library('MatrixGenerics')"

RUN R -e  "BiocManager::install('SummarizedExperiment')"
RUN R -e "library('SummarizedExperiment')"

RUN R -e  "BiocManager::install('BiocFileCache')"
RUN R -e "library('BiocFileCache')"


RUN R -e  "BiocManager::install('SingleCellExperiment')"
RUN R -e "library('SingleCellExperiment')"


# RUN apt-get update && apt-get install libglib2.0-0t64=2.80.2-2  libglib2.0-bin=2.80.2-2  -y --allow-downgrades
RUN apt-get update && apt-get install libglib2.0-dev-bin  -y --allow-downgrades


RUN apt-get update && apt-get install libmagick++-dev  -y --allow-downgrades

# RUN R -e "if ('SpatialExperiment' %in% rownames(installed.packages())) remove.packages('SpatialExperiment')"

RUN R -e "BiocManager::install('SpatialExperiment', ask=FALSE)"
RUN R -e "library('SpatialExperiment')"


RUN R -e "install.packages('remotes')"
RUN R -e "remotes::install_version('Matrix', version = '1.6.0')"

RUN R -e "BiocManager::install('AnnotationDbi')"
RUN R -e "library('AnnotationDbi')"

RUN R -e "BiocManager::install('edgeR')"
RUN R -e "library('edgeR')"

RUN R -e "BiocManager::install('genefilter')"
RUN R -e "library('genefilter')"

RUN R -e "BiocManager::install('BiocParallel')"
RUN R -e "library('BiocParallel')"

RUN R -e "remotes::install_version('Matrix')"

RUN R -e "BiocManager::install('zinbwave')"
RUN R -e "library('zinbwave')"

#s  cran dep
RUN R -e "BiocManager::install('metapod')"
RUN R -e "library('metapod')"

RUN R -e "BiocManager::install('BiocNeighbors')"
RUN R -e "library('BiocNeighbors')"

RUN R -e "BiocManager::install('bluster')"
RUN R -e "library('bluster')"

RUN R -e "BiocManager::install('beachmat')"
RUN R -e "library('beachmat')"

RUN R -e "BiocManager::install('ScaledMatrix')"
RUN R -e "library('ScaledMatrix')"

RUN R -e "BiocManager::install('DelayedMatrixStats')"
RUN R -e "library('DelayedMatrixStats')"

RUN R -e "BiocManager::install('BiocSingular')"
RUN R -e "library('BiocSingular')"

RUN R -e "BiocManager::install('limma')"
RUN R -e "library('limma')"

RUN R -e "BiocManager::install('scran')"
RUN R -e "library('scran')"

RUN apt-get update && apt-get install cmake -y
RUN R -e "BiocManager::install('SpatialExperiment')"
RUN R -e "library('SpatialExperiment')"

RUN R -e "install.packages('SpatialDDLS')"
RUN R -e "library('SpatialDDLS')"

RUN R -e  "install.packages('R.utils')"
RUN R -e "library('R.utils')"


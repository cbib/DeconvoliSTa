# base image 
FROM pytorch/pytorch:2.1.0-cuda11.8-cudnn8-runtime

# install os packages
RUN apt-get update && apt-get install -y git

# install miniconda3 - https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html
RUN conda create -y -n cell2loc_env python=3.10 pip

# create conda environment yaml file
RUN /opt/conda/bin/conda run -n cell2loc_env pip install --no-cache-dir \
        numpy==1.21.4 \
        pandas==1.3.4 \
        scipy==1.7.2 \
        scikit-learn==1.0.1 \
        anndata==0.8.0 \
        scanpy==1.8.2 \
        pytorch-lightning==1.3.8 \
        scvi-tools==0.14.3 \
        cell2location==0.6a0
    /opt/conda/bin/conda clean -afy
ENV PATH /opt/conda/bin:$PATH


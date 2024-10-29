#  DeconvolisSTa **Deconvol**ut**i**on of **S**patial **T**ranscriptomics d**A**ta

## Introduction

This document provides documentation for the Snakemake version of Spotless, a spatial deconvolution pipeline. The pipeline was developed during my internship at CBIB. Deconvolution can be performed using one or more of the following methods: Cell2location, RCTD, NNLS, SpatialDWLS, DDLS, Seurat, MusiC, and Dirichlet (random).

## Required environment installation

Here are needed utilities softwares to run the pipeline.

### Snakemake

The pipeline is implemented in Snakemake which can be installed via the tutorial available in this link [Setting up Conda and Snakemake](https://gist.github.com/RomainFeron/da9df092656dd799885b612fedc9eccd).

### Singularity

Almost all available methods for deconvolution are containerized and run in docker images. However, Snakemake is compatible only with Singularity. Snakemake converts the docker image to a singularity image '.sif' before using it. To install Singularity, this [tutorial](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) from official documentation is useful.

### System requirements

The CPU used to execute the pipeline is an `Intel(R) Xeon(R) CPU E5-1607 0 @ 3.00GHz` model with 4 nodes and `40 Gi` RAM. And the GPU used is `NVIDIA Quadro P6000 24GB PCIe 3.0` with Cuda driver version `535.171.04` and Cuda `12.2`.

With a single cell file of 3.3 Gi, this machine couldn't run the pipeline because of an out of memory. So, think about optimizing single cell data file size or dividing it into chunks as we did for ours.

Here are some execution times for synthetic data.

| | 100 spots | 1000 spots | 10000 spots |
|-|-----------|------------|-------------|
| Cell2location | 2h 16m | 3h 55m | 3h 37m |
| RCTD | 7m | 41m | 1h |

## Pipeline running

Here is how to run the pipeline.

### Running parameters

The pipeline inputs are as following. First, `sc_input` the singleCell reference file and `sp_input` the spatial transcriptomic file. The two files should be in RDS formats. The singleCell file should have a column for annotation in its `meta.data` attribute. These annotations will be used as cell types for the deconvolution.  Its is recommended to verify existence of this column before running the pipeline. If the single cell file is too big (> 3 Go), it is recommended to reduce its size be removing some unnecessary data for deconvolution, for example deleting unused `annotation_level` in case of multiple `annotation_levels`. Or, chunking the file into equal cell proportions subsamples.

Spatial and single cell files should have the same gene barcoding type, if it is not the case, a mapping should be done between gene names for the two files to homogenize gene names.

Here's an example of the command line to run the pipeline:

```bash
snakemake  -s   main.smk -c8 --config \
    mode="run_dataset"  methods=cell2location,rctd \
    sc_input="test_sc_data.rds" sp_input="test_sp_data.rds" \
    output="res" use_gpu="true" skip_metrics="true" \
    annot="subclass" map_genes="false" load_model="true" \
    model_path="mod" --use-singularity \
    --singularity-args '\--nv'
```


1. The number of cores to use for Snakemake is specified with -c8, the parameters for Snakemake are specified in a dictionary named config.

2. methods is the list of methods to run, pay attention to separate them with a comma without spaces.

3. output is the directory where output files will be.

4. use_gpu is a parameter to use a GPU running. By default, the CPU is used.

5. if skip_metrics is set to true, benchmarking metrics are not performed, for instance correlation, RMSE , accuracy, balanced accuracy, sensitivity. If you do want to measure them you should have a data frame of deconvolution ground truth in the $relative_spot_composition attribute in the spatial data input sp_input. The default value of this parameter is false.

6. annot is the name of annotation column in the single cell input metadata. Pay attention to verify that it is a valid column name before running. In addition, this attribute should be a char vector not a factor. The default value of this parameter is subclass.

7. map_genes should be set to true if the gene names in the single cell and spatial data inputs are not in the same gene symbols format. The default value of this parameter is false.
8. ---use-singularity and ---singularity-args '---nv' to enable singularity use and GPU access for Snakemake.

9. When load_model is true, the cell2location model doesn't do the build stage in the pipeline, instead it is loaded from model_path. When having multiple spatial samples associated with the same single cell reference dataset, this feature allows to do the build of cell2location model once and do the predictions for all spatial samples without rebuilding the model each time.

# Synthetic data generation with Synthspot
## With the pipeline, you can generate synthetic spatial data. Taking a single cell data input, Synthspot generate different synthetic spatial profiles. The generated datasets types include a variety of synthetic and real-world data configurations designed to capture different spatial and cell type distributions.

## The listing below shows an example of line command to generate an artificial dataset from golden standard dataset provided by spotless.

```bash
  snakemake -s  main.smk -c12 --config mode="generate_data" \
      sc_input="standards/reference/gold_standard_1.rds" \
      dataset_type="aud" reps=1  output="synthetic_data_sm" \
      region_var="celltype_coarse" --use-singularity

```
1. sc_input is the single cell input file used to generate synthetic data.
2. dataset_type is the dataset profile to generate. This parameter is a one or a comma separated list of types mapped in the listing below.
3. rep is the number of replicates generated for the generated dataset types selected. The output will be $N_{types}\times rep$ files indexed with suffix {type}_{rep}.rds.
4. output is the output directory for generated files.
5. region_var column with regional metadata in sc_input@meta.data, if any (for "real" dataset types).


# Spatial Transcriptomics Visualizations
### The pipeline include an interactive visualization tool for deconvolution results. It is an independent part of the pipeline that shows deconvolution results with proportions for each spot displayed with the actual Visium spatial image, and different deconvolution methods results can be visualized. In addition, it displays spots with clustering. The tool include also raw visualized data. A demo of this tool can be seen   <a href="https://drive.google.com/uc?export=download&id=1eXaHzJOT6B9YIPYDvQtTKoAs0kv8eoiX" download target="_blank" rel="noopener noreferrer">here</a>


```bash
  snakemake -s  main.smk -c8 --config mode="generate_vis" \
  sp_input="UKF243_T_ST_1_raw.rds" output="vis_output" \
  norm_weights_filepaths="props_rctd.tsv,props_cell2location.tsv" \
  st_coords_filepath="tissue_positions_list_243.csv" \
  data_clustered="seurat_metadata.csv" image_path="tissue_hires.png" \
  scale_factor='0.24414062' deconv_methods=rctd,cell2location
```
1. sp_input is the spatial data file used in deconvolution, it is used only to name the output HTML file.
2. output is the output directory.
3. norm_weights_filepaths is a list of deconvolution results files when using multiple deconvolution methods. The filenames should be comma separated without spaces in between. In addition, files should be tabulation separated values (TSV) files.
4. st_coords_filepath is a CSV file having the spots as index and their correspondent pixel positions in the Visium image, in pxl_col_in_fullres and pxl_row_in_fullres columns. The file should have six columns as shown in the example below. The preprocessing of data add columns names in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres.

| Spot | in_tissue | array_row | array_col | pxl_row_in_fullres | pxl_col_in_fullres |
|------|-----------|-----------|-----------|---------------------|---------------------|
| ACGCCTGACACGCGCT-1 | 0 | 0 | 0 | 721 | 1375 |
| TACCGATCCAACACTT-1 | 0 | 1 | 1 | 796 | 1418 |
| ATTAAAGCGGACGAGC-1 | 0 | 0 | 2 | 721 | 1461 |

5. data_clustered is another CSV file which associate every spot to its cluster using a column named BayesSpace.
6. image_path is he path to Visium image of the sample.
7. scale_factor is the scaling factor to use in order to match pixels coordinates with the Visium image.
8. deconv_methods is a list of deconvolution methods used in the norm_weights_filepaths files.

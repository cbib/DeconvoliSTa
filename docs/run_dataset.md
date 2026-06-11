# Mode: `run_dataset`

Deconvolve a spatial transcriptomics dataset. For each method in `methods`, the pipeline writes a
proportions table `proportions_<method>_<sample>.tsv` (rows = spots, columns = cell types, each row
sums to 1). Optionally it also builds the interactive visualization (`do_visu`) and/or evaluation
metrics (`skip_metrics=false`, needs a ground-truth dataset).

This is the main entry point. See `SETUP.md` for the full getting-started guide.

## Required configuration (`--config ...` or `config.yaml`)

| Key       | Example                       | Notes |
|-----------|-------------------------------|-------|
| `mode`    | `"run_dataset"`               | |
| `methods` | `"rctd,nnls,cell2location"`   | any of `rctd, cell2location, nnls, spatialdwls, ddls, dirichlet` (comma-separated) |
| `sc_input`| `"reference.rds"`             | single-cell reference (Seurat `.rds`); used by every method **except** `dirichlet` |
| `sp_input`| `"sample.rds"`                | spatial Visium sample (Seurat `.rds`) |
| `annot`   | `"cell_type"`                 | name of the metadata column in the reference holding the cell-type labels |

## Optional configuration

| Key            | Default  | Notes |
|----------------|----------|-------|
| `output`       | `.`      | output directory |
| `map_genes`    | `"false"`| `"true"` converts gene symbols ↔ Ensembl when the reference and spatial use different gene IDs (needs the `org.Hs.eg.db` images for nnls / spatialdwls — see SETUP.md) |
| `do_visu`      | `"false"`| `"true"` → after deconvolution, extract the tissue image/coords from the spatial object, cluster the spots (BayesSpace + Seurat) and build the interactive HTML (needs `sp_bayesspace.sif` + `visu.sif`) |
| `bayes_q`      | `"auto"` | number of BayesSpace spatial domains: `"auto"` (qTune, slower) or an integer, e.g. `"9"` |
| `seurat_res`   | `"0.8"`  | Seurat clustering resolution (higher → more clusters) |
| `n_largest_cell_types` | `"5"` | cell types shown per pie in the visualization |
| `skip_metrics` | `"false"`| `"true"` to skip the evaluation metrics |
| `sif_dir`      | `"sif"`  | directory searched for local `.sif` images (else the public image is pulled) |

## Examples

Deconvolution only (public images pulled automatically, nothing to install):
```bash
snakemake -s main.smk -c4 --use-singularity \
  --config mode="run_dataset" methods="rctd,nnls,spatialdwls" \
  sc_input="reference.rds" sp_input="sample.rds" \
  output="results" annot="cell_type"
```

Full pipeline, all the way to the interactive visualization:
```bash
snakemake -s main.smk -c4 --use-singularity \
  --config mode="run_dataset" methods="rctd,nnls,spatialdwls,cell2location" \
  do_visu="true" sif_dir="$SIF_DIR" \
  sc_input="reference.rds" sp_input="sample.rds" \
  output="results" annot="cell_type"
```

## Output
```
results/
├── proportions_<method>_<sample>.tsv      # one per method
├── metrics/ ...                            # if skip_metrics=false
└── visualization/                          # if do_visu=true
    ├── <sample>_visualization.html         # the interactive report
    ├── assets/ (tissue_image.png, tissue_positions.csv, scale_factor.txt)
    ├── clustering_bayesspace.csv
    └── clustering_seurat.csv
```

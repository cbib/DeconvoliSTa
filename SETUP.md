# DeconvoliSTa — setup & run

Reproducible benchmarking pipeline for spatial transcriptomics deconvolution
(Snakemake + Singularity). It **runs locally on a single machine** — no cluster
required — and can *optionally* be submitted to a SLURM HPC cluster for large or
GPU jobs.

## 1. Prerequisites

- **Snakemake** in a conda environment (e.g. `conda create -n snakemake -c bioconda -c conda-forge snakemake`).
- **Singularity / Apptainer** (runs the per-method container images, including `docker://` ones).
- That's all for a local run. *Optional:* a **SLURM** cluster (a CPU partition, and a GPU
  partition for the `cell2location` / `ddls` GPU runs) if you want to submit batch jobs.

## 2. Configuration (no hardcoded paths)

All site-specific paths live in a single file. Copy the template and edit it once:

```bash
cp config.env.example config.env
nano config.env        # set repo dir, output dir, SIF dir, conda, partitions, GPU resource
```

`config.env` is gitignored — your personal paths never reach the repo.

## 3. Container images

Most methods pull their image automatically via `docker://` (Snakemake converts to a SIF on first use).
A few images were rebuilt locally to fix missing dependencies / CUDA compatibility; their build
recipes are versioned next to each method:

| Image (in `$SIF_DIR`)         | Built from                                                  |
|-------------------------------|-------------------------------------------------------------|
| `sp_nnls_cbib.sif`            | `subworkflows/deconvolution/nnls/nnls.def`                  |
| `sp_spatialdwls_cbib.sif`     | `subworkflows/deconvolution/spatialdwls/spatialdwls.def`    |
| `sp_cell2location_cu128.sif`  | `subworkflows/deconvolution/cell2location/cell2location_cu128.def` |
| `sp_rctd_latest.sif`          | provided image (place a copy in `$SIF_DIR`)                 |

Build one (needs internet + `--fakeroot`, e.g. on a build host):

```bash
singularity build --fakeroot $SIF_DIR/sp_nnls_cbib.sif subworkflows/deconvolution/nnls/nnls.def
```

Make sure every `.sif` listed above sits in the `$SIF_DIR` you set in `config.env`.

## 4. Run

### Quick local run (no cluster)

Run any method directly with Snakemake; Singularity pulls the image on first use.
The bundled unit-test data runs in seconds and needs no configuration:

```bash
snakemake -s main.smk -c4 --use-singularity \
  --config mode="run_dataset" methods="rctd" \
  sc_input="unit-test/test_sc_data.rds" \
  sp_input="unit-test/test_sp_data.rds" \
  output="res" annot="subclass"
```

Swap `methods=` for any of `rctd, nnls, spatialdwls, dirichlet, cell2location, ddls`
(comma-separated to run several). Results land in `res/proportions_<method>_*.tsv`.
For your own data, point `sc_input` / `sp_input` at your `.rds` files and set `annot`
to the cell-type column of your single-cell reference.

### On a SLURM cluster (optional)

Use the wrapper — it reads `config.env` and submits the right SLURM job:

```bash
./launch.sh run_deconvolista.sh          # CPU batch, all methods, on the test dataset
./launch.sh run_cell2location_gpu.sh     # cell2location on GPU
./launch.sh run_ddls_gpu.sh              # DDLS on GPU
```

Override inputs/methods on the fly (forwarded to the job via `--export=ALL`):

```bash
METHODS=rctd,nnls,dirichlet ./launch.sh run_deconvolista.sh
SC_INPUT=my_ref.rds SP_INPUT=my_sample.rds ANNOT=cell_type ./launch.sh run_deconvolista.sh
```

Results are written to `$OUTPUT_DIR` (one `proportions_<method>_*.tsv` per method, plus
comparison metrics when `skip_metrics=false`).

## 5. Interactive visualization

See `docs/generate_vis.md` — builds an interactive HTML overlaying deconvolution proportions on
the Visium image with spot clustering.

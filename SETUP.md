# DeconvoliSTa — setup & run

Reproducible benchmarking pipeline for spatial transcriptomics deconvolution
(Snakemake + Singularity). It **runs locally on a single machine**, no cluster required and can *optionally* be submitted to a SLURM HPC cluster for large or GPU jobs.

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

**You do not need to build anything for the deconvolution.** Each method runs in a container that
Snakemake pulls automatically from a public registry on first use (`--use-singularity`).

Build a local `.sif` only for the cases below, and drop it in the `$SIF_DIR` from `config.env`
(the pipeline uses a local image when it is present, otherwise it falls back to the public one):

| Build only if you use… | Image(s) | Build from |
|------------------------|----------|------------|
| `map_genes="true"` (convert gene symbols ↔ Ensembl, for nnls / spatialdwls) | `sp_nnls_cbib.sif`, `sp_spatialdwls_cbib.sif` | `subworkflows/deconvolution/nnls/nnls.def`, `.../spatialdwls/spatialdwls.def` |
| `cell2location` on a recent GPU (CUDA 12.8 / 570-series driver) | `sp_cell2location_cu128.sif` | `subworkflows/deconvolution/cell2location/cell2location_cu128.def` |
| `do_visu="true"` (BayesSpace clustering + interactive HTML) | `sp_bayesspace.sif`, `visu.sif` | `subworkflows/clustering/bayesspace.def`, `subworkflows/visualization/visu.def` |

```bash
# example (needs internet + --fakeroot, e.g. on a build host)
singularity build --fakeroot $SIF_DIR/sp_bayesspace.sif subworkflows/clustering/bayesspace.def
```

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

`launch.sh` reads `config.env` and submits the right SLURM job (correct partition, GPU, logs).

**The easy way — `run_all.sh` (you choose the methods, it does the rest):**
```bash
SC_INPUT=my_ref.rds SP_INPUT=my_sample.rds ANNOT=cell_type ./run_all.sh
```
Pick any set of methods with `METHODS=` and `run_all` **routes each one to the right hardware**
(CPU vs GPU), runs them **in parallel**, then submits a **visualization job that waits for them all**
(SLURM `--dependency=afterok`) and builds the HTML over exactly those methods:
```bash
METHODS=cell2location            ... ./run_all.sh   # just one method
METHODS=rctd,cell2location       ... ./run_all.sh   # any subset
DO_VISU=false METHODS=rctd       ... ./run_all.sh   # proportions only, no visualization
```
Knobs: `METHODS` (default `rctd,nnls,spatialdwls,cell2location,ddls`), `DO_VISU` (default `true`),
`BAYES_Q`. Track with `squeue -u $USER`.

**Or submit jobs individually** (a single SLURM job runs on a single partition, so CPU and GPU
methods are separate scripts):
```bash
DO_VISU=true METHODS=rctd,nnls,spatialdwls ./launch.sh run_deconvolista.sh  # CPU methods (+ visu)
./launch.sh run_cell2location_gpu.sh                                        # cell2location on GPU
./launch.sh run_ddls_gpu.sh                                                 # ddls on GPU
```

Env vars are forwarded to the job (`--export=ALL`), so any config option can be overridden on the
fly (`METHODS`, `SC_INPUT`, `SP_INPUT`, `ANNOT`, `DO_VISU`, `BAYES_Q`, …).

Results are written to `$OUTPUT_DIR` (one `proportions_<method>_*.tsv` per method, plus the
`visualization/` folder when the visualization runs, and metrics when `skip_metrics=false`).

## 5. Interactive visualization

### Automatic — one command does everything (`do_visu=true`)

Add `do_visu="true"` to a `run_dataset` run and the pipeline goes all the way to the interactive
HTML **on its own**, with no manually-prepared files:

```bash
snakemake -s main.smk -c4 --use-singularity \
  --config mode="run_dataset" methods="rctd,nnls,spatialdwls,cell2location" \
  do_visu="true" \
  sc_input="reference.rds" sp_input="sample.rds" \
  output="results" annot="cell_type" \
  sif_dir="$SIF_DIR"
```

After the deconvolution, it automatically:
1. **extracts** the tissue image, spot coordinates and scale factor from the spatial Seurat
   object (`@images`) — nothing to prepare by hand;
2. computes **two clusterings** of the spots: BayesSpace (spatial domains) and Seurat
   (transcriptomic) — the visualization shows a toggle to switch between them;
3. builds the **interactive HTML**.

Output: `results/visualization/<sample>_visualization.html` (plus the extracted assets and the two
clustering CSVs next to it).

**Requirements (visualization only).** The deconvolution itself needs no custom image (public
images are pulled automatically). The visualization step uses **two images that are not on a public
registry** — build them once and put them in `$SIF_DIR`:
```bash
singularity build --fakeroot $SIF_DIR/sp_bayesspace.sif subworkflows/clustering/bayesspace.def
singularity build --fakeroot $SIF_DIR/visu.sif           subworkflows/visualization/visu.def
```
If the spatial object has **no image** in `@images`, the extraction fails on purpose — run with
`do_visu="false"` (the default) to get the proportions only.

**Optional knobs** (config):
- `bayes_q` — number of BayesSpace spatial domains; `"auto"` (default, runs qTune, slower) or an
  integer, e.g. `bayes_q="9"`.
- `seurat_res` — Seurat clustering resolution (default `0.8`).
- `n_largest_cell_types` — cell types shown per pie (default `5`).

On the cluster, `run_all.sh` / `launch.sh` forward `do_visu` like any other override:
```bash
DO_VISU=true SC_INPUT=ref.rds SP_INPUT=sample.rds ANNOT=cell_type ./launch.sh run_deconvolista.sh
```

### Manual (advanced)

`mode="generate_vis"` builds the HTML from already-computed proportions + an image + coordinates +
a clustering you provide yourself. See `docs/generate_vis.md`. The automatic mode above is the
recommended path.

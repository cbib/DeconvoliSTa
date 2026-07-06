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

**You do not need to build anything.** Every step — deconvolution, clustering and the
visualization — runs in a container that Snakemake pulls automatically on first use
(`--use-singularity`). Most images are on the CBiB registry (`ghcr.io/cbib/*`), a few helpers on
Docker Hub. The cbib R images already bundle `org.Hs.eg.db` (so `map_genes="true"` works) and the
cell2location image ships a torch CUDA 12.8 build for recent GPUs — nothing to build for those
cases. See `docs/data_and_resources.md` for the full image list.

> The `ghcr.io/cbib` packages must be **public** for an anonymous pull (otherwise run
> `singularity remote login` / `docker login ghcr.io` once). A local-`.sif` override via `$SIF_DIR`
> is planned but not wired in yet; images are pulled from the registry for now. The `.def` /
> `Dockerfile.cbib` recipes are kept for rebuilding an image.

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

### With the `./deconvolista` launcher (SLURM or local)

A single launcher, **`./deconvolista`**, reads `config.env` and runs everything for you: it
picks any set of methods, **routes each to the right hardware** (CPU vs GPU), runs the GPU methods
**in parallel**, then runs the **visualization after they all succeed** (on SLURM via a
`--dependency=afterok` job) and builds the HTML over exactly those methods.

On a cluster it submits SLURM jobs; if `sbatch` is not available (or with `--local`) it runs here
directly in one `snakemake` call. **Before running anything**, it prints a full execution plan
(inputs, CPU/GPU routing, SLURM resources, the `snakemake` command) and asks for confirmation —
preview it with `--dry-run`, or skip the prompt in scripts with `--yes`.

**Interactive (just run it):**
```bash
./deconvolista                 # menu: pick methods, visualization, q; it shows the equivalent command
```

**Non-interactive (scriptable):**
```bash
./deconvolista --methods cell2location                  # just one method
./deconvolista --methods rctd,cell2location --visu      # any subset, with the visualization
./deconvolista --methods rctd --no-visu                 # proportions only
./deconvolista --local --methods rctd                   # run here, no SLURM
./deconvolista --dry-run --methods nnls,ddls            # print the execution plan, run nothing
./deconvolista --yes --methods rctd,nnls                # skip the confirmation prompt (scripts)
./deconvolista --help
```
Inputs come from `config.env` (`SC_INPUT` / `SP_INPUT` / `ANNOT`); override per run with
`--sc` / `--sp` / `--annot`. Other knobs: `--bayes-q`, and `DEFAULT_METHODS` / `DO_VISU` in
`config.env`. Track with `squeue -u $USER`.

> Methods are **auto-discovered**: any `subworkflows/deconvolution/<m>/run_<m>.smk` is a valid
> `--methods <m>`, and its `method.env` (`USE_GPU=true|false`) decides CPU vs GPU. Adding a method
> needs no edit to the launcher. Resources per job type are set in `config.env`
> (`MEM_CPU` / `MEM_GPU` / `MEM_VISU` / `TIME_LIMIT`).

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

**Requirements (visualization only).** Nothing to build — the two visualization images
(`ghcr.io/cbib/deconvolista-bayesspace` for BayesSpace + Seurat clustering, and
`ghcr.io/cbib/deconvolista-visu` for the HTML) are pulled automatically like the others.
If the spatial object has **no image** in `@images`, the extraction fails on purpose — run with
`do_visu="false"` (the default) to get the proportions only.

**Optional knobs** (config):
- `bayes_q` — number of BayesSpace spatial domains; `"auto"` (default, runs qTune, slower) or an
  integer, e.g. `bayes_q="9"`.
- `seurat_res` — Seurat clustering resolution (default `0.8`).
- `n_largest_cell_types` — cell types shown per pie (default `5`).

On the cluster, `./deconvolista` builds the visualization by default (`--visu`); use `--no-visu`
to skip it:
```bash
./deconvolista --methods rctd,nnls,spatialdwls --visu
```

### Manual (advanced)

`mode="generate_vis"` builds the HTML from already-computed proportions + an image + coordinates +
a clustering you provide yourself. See `docs/generate_vis.md`. The automatic mode above is the
recommended path.

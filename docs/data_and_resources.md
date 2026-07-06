# Data and Resources

Overview of data and compute resources needed to run and test the DeconvoliSTa pipeline.

---

## Test data (unit tests)

Minimal dataset for validating that the pipeline runs correctly. Suitable for CI and quick local checks.

| File | Location | Description |
|------|----------|-------------|
| `test_sc_data.rds` | `unit-test/` (repo) | Single-cell reference — 26,252 genes × 1,404 cells. Derived from Allen cortex dataset, genes with < 10 counts filtered out. Annotation column: `subclass` |
| `test_sp_data.rds` | `unit-test/` (repo) | Spatial data — 26,252 genes × 16 spots. Synthetically generated from `test_sc_data` with `synthspot` (`artificial_uniform_distinct`, seed 10) |

**Example run with test data:**
```bash
snakemake -s main.smk --config mode="run_dataset" \
  methods="rctd" \
  sc_input="unit-test/test_sc_data.rds" \
  sp_input="unit-test/test_sp_data.rds" \
  output="res" annot="subclass" \
  --use-singularity
```

---

## Gold standard data

Real spatial transcriptomics dataset with known ground truth, used for scientific benchmarking.

| File | Location | Description |
|------|----------|-------------|
| `gold_standard_1.rds` | `standards/reference/` (repo) | Single-cell reference |
| `Eng2019_cortex_svz_fov*.rds` | `standards/gold_standard_1/` (repo) | 7 spatial samples (seqFISH, cortex/SVZ) |
| `Eng2019_cortex_svz_fov*.jpeg` | `standards/gold_standard_1/` (repo) | Tissue images for visualisation |

**Example run with gold standard:**
```bash
snakemake -s main.smk --config mode="run_dataset" \
  methods="rctd" \
  sc_input="standards/reference/gold_standard_1.rds" \
  sp_input="standards/gold_standard_1/Eng2019_cortex_svz_fov5.rds" \
  output="res" annot="celltype" \
  --use-singularity
```

> Larger real datasets (e.g. a single-cell reference + Visium sample) are not shipped with the
> repo. Point `sc_input`/`sp_input` at your own `.rds` files; set `annot` to the metadata column
> holding the cell-type labels of your reference.

---

## Container images

Every step runs in a Singularity/Apptainer container that Snakemake **pulls automatically** on
first use (`--use-singularity`) — nothing to build. Most images live on the CBiB GitHub Container
Registry (`ghcr.io/cbib/*`); a few helpers are on Docker Hub.

| Step | Image |
|------|-------|
| RCTD | `docker://ghcr.io/cbib/sp_rctd:latest` |
| NNLS | `docker://ghcr.io/cbib/sp_nnls:latest` |
| spatialDWLS | `docker://ghcr.io/cbib/sp_spatialdwls:latest` |
| cell2location | `docker://ghcr.io/cbib/sp_cell2location:latest` (GPU, torch CUDA 12.8) |
| DDLS | `docker://ghcr.io/cbib/sp_ddls:latest` (GPU) |
| Dirichlet | `docker://csangara/seuratdisk:latest` |
| BayesSpace clustering (`do_visu`) | `docker://ghcr.io/cbib/deconvolista-bayesspace:latest` |
| Visualization (`do_visu`) | `docker://ghcr.io/cbib/deconvolista-visu:latest` |
| RDS ↔ H5AD conversion | `docker://csangara/seuratdisk:latest` |
| Evaluation (`skip_metrics=false`) | `docker://csangara/sp_eval:latest` |
| Synthetic data (`generate_data`) | `docker://csangara/synthspot:latest` |

The `ghcr.io/cbib/*` R images bundle `org.Hs.eg.db`, so `map_genes="true"` works out of the box;
the cell2location image ships a torch CUDA 12.8 build for recent GPUs (570-series driver).

> The `ghcr.io/cbib` packages must be **public** for an anonymous pull to succeed (otherwise run
> `singularity remote login` / `docker login ghcr.io` first).
>
> The `.def` / `Dockerfile.cbib` recipes next to each method are kept for rebuilding an image. A
> local-`.sif` override (`SIF_DIR` in `config.env`) is planned but not wired in yet — images are
> pulled from the registry for now.

---

## Compute resources (rough estimates)

Indicative figures for sizing SLURM jobs on an HPC cluster. Adjust to your data and hardware.

| Method | Job type | Cores | RAM (est.) | Runtime (est.) |
|--------|----------|-------|-----------|----------------|
| RCTD | CPU | 12 | 32 GB | 7–60 min |
| cell2location | GPU | 8 | 16 GB | 2–4 h |
| DDLS | GPU | 12 | 32 GB | 4–12 h |
| NNLS | CPU | 8 | 8 GB | < 10 min |
| spatialDWLS | CPU | 12 | 16 GB | 15–30 min |
| Dirichlet | CPU | 8 | 8 GB | < 10 min |

Runtimes measured on real data (~1000 spots, ~10,000 genes). Unit-test data runs in seconds.

GPU methods (cell2location, DDLS) need an NVIDIA GPU and the container's CUDA build to match the
host driver. If your cluster uses H100 MIG slices, request one via `GPU_GRES` in `config.env`
(e.g. `gpu:nvidia_h100_nvl_1g.24gb`); a 24 GB slice is enough for the demo data.

---

## Known limitations / TODO

- `seurat` is listed in `main.smk`'s supported methods but has no `run_seurat.smk` — it will fail if selected.
- `ICA` and `CytoSPACE` are not integrated yet (planned).

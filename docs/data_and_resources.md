# Data and Resources

Overview of data and compute resources needed to run and test the DeconvoliSTa pipeline.

---

## Test data (unit tests)

Minimal dataset for validating that the pipeline runs correctly. Suitable for CI and quick local checks.

| File | Location | Description |
|------|----------|-------------|
| `test_sc_data.rds` | `unit-test/` (repo) | Single-cell reference — 26,252 genes × 1,404 cells. Derived from Allen cortex dataset, genes with < 10 counts filtered out. Annotation column: `subclass` |
| `test_sp_data.rds` | `unit-test/` (repo) | Spatial data — 26,252 genes × 16 spots. Synthetically generated from `test_sc_data` with `synthspot` (`artificial_uniform_distinct`, seed 10) |

Also available on Apollo at `/mnt/cbib/RetinRNA/spatial/` (`.rds` and `.h5ad` formats).

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

---

## Singularity images

Pre-built container images required to run each method. Stored on Apollo at `/mnt/cbib/RetinRNA/spatial/`.

| Method | Image on Apollo | Docker Hub fallback |
|--------|----------------|---------------------|
| RCTD | `sp_rctd_latest.sif` ✓ | `docker://csangara/sp_rctd:latest` |
| cell2location | `sp_cell2location_h100.sif` ✓ | `docker://csangara/sp_cell2location:latest` |
| spatialDWLS | — | `docker://csangara/sp_spatialdwls:latest` |
| NNLS | — | `docker://csangara/sp_nnls:latest` |
| DDLS | — | `docker://csangara/sp_ddls:latest` |
| Dirichlet | — | `docker://csangara/seuratdisk:latest` *(image de conversion utilisée en attendant une image dédiée)* |
| Seurat | — | *(no .smk file yet — not runnable)* |
| Evaluation | — | `docker://csangara/sp_eval:latest` |
| RDS↔H5AD conversion | — | `docker://csangara/seuratdisk:latest` |

> **Note:** Methods without a local `.sif` file will pull from Docker Hub on first run, which is slow and requires internet access. Building local `.sif` files for all methods is recommended for production use on Apollo.

---

## Compute resources (Apollo cluster)

| Method | Partition | Cores | RAM (est.) | Runtime (est.) |
|--------|-----------|-------|-----------|----------------|
| RCTD | `compute` | 12 | 32 GB | 7–60 min |
| cell2location | `gpu` | 8 | 16 GB | 2–4 h |
| DDLS | `compute` | 12 | 32 GB | 4–12 h |
| NNLS | `compute` | 8 | 8 GB | < 10 min |
| spatialDWLS | `compute` | 12 | 16 GB | 15–30 min |
| Dirichlet | `compute` | 8 | 8 GB | < 10 min |

Runtimes measured on real data (~1000 spots, ~10,000 genes). Unit test data runs in seconds.

For cell2location, available GPU MIG types on Apollo:
- `nvidia_h100_nvl_4g.47gb` — recommended
- `nvidia_h100_nvl_1g.24gb` — minimum viable

---

## Known limitations / TODO

- `seurat` is listed as a supported method in `main.smk` but has no `run_seurat.smk` — will fail if selected
- `use_gpu` is hardcoded to `"false"` in `run_cell2location.smk` — GPU support currently disabled
- No local `.sif` for spatialDWLS, NNLS, DDLS, Dirichlet, Seurat — Docker Hub pull required

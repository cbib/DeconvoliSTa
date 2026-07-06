## Method parameters

Methods integrated in the pipeline:
- [cell2location](#cell2location)
- [spatialDWLS](#spatialdwls)
- RCTD — no adjustable parameters
- NNLS — no adjustable parameters
- Dirichlet — no adjustable parameters (random baseline / negative control)
- DDLS — epochs / batch size set in `run_ddls.smk` (see below)

> `seurat` has a folder but **no `run_seurat.smk`** — it is not wired into the pipeline yet and will
> fail if selected.

### Usage

Tunable parameters are provided in `config.yaml` under the `deconv_args` key, with the method name
**in lowercase**:

```yaml
# config.yaml
deconv_args:
  cell2location:
    build: "-e 250"
    fit: "-n 8 -d 200"
  spatialdwls: "--n_topmarkers 50"
```

Then run:
```bash
snakemake -s main.smk --use-singularity --config mode="run_dataset" \
  methods=cell2location,spatialdwls \
  sc_input="path/to/sc.rds" sp_input="path/to/sp.rds" output="results" annot="cell_type"
```

**DDLS** does not read `deconv_args`; its `--epochs` and `--batch_size` are set directly in
`subworkflows/deconvolution/ddls/run_ddls.smk`.

---

#### cell2location

Model building (`deconv_args.cell2location.build`):
- `-s`: metadata column in the single-cell reference containing sample / batch information (default: None)
- `-t`: metadata column containing multiplicative technical effects, e.g. platform effects (default: None)
- `-e`: number of epochs to train the model (default: 250)
- `-p`: number of samples to take from the posterior distribution (default: 1000)

Model fitting (`deconv_args.cell2location.fit`):
- `-n`: estimated number of cells per spot (default: 8)
- `-d`: within-experiment variation in RNA detection sensitivity (default: 200)
- `-e`: number of epochs to fit the model (default: 30000)
- `-p`: number of samples to take from the posterior distribution (default: 1000)

#### spatialDWLS
- `--n_topmarkers`: number of top marker genes per cell type to use (default: 100)
- `--nn.dims`: number of PCs to use for the nearest-network calculation (default: 10)
- `--nn.k`: number of neighbors for the nearest-network calculation (default: 4)
- `--cluster.res`: Leiden cluster resolution (default: 0.4)
- `--cluster.n_iter`: number of iterations during Leiden clustering (default: 1000)

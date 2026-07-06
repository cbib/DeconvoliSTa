# Adding a deconvolution method

Methods are **auto-discovered** — adding one needs **no edit to `main.smk` or `./deconvolista`**.
`main.smk` treats `<method>` as valid iff `subworkflows/deconvolution/<method>/run_<method>.smk`
exists, and `./deconvolista` picks it up automatically and routes it to CPU or GPU from its
`method.env`.

## 1. Create the method folder

`subworkflows/deconvolution/<method>/` containing:

| File | Required | Purpose |
|------|----------|---------|
| `run_<method>.smk` | ✅ | Snakemake rule that produces the proportions table |
| `method.env` | ✅ | one line — `USE_GPU=true` or `USE_GPU=false` (read by `./deconvolista` to route the job to the GPU or CPU partition) |
| `script.R` / `script.py` | ✅ | your actual algorithm |
| `Dockerfile.cbib` / `<method>.def` | optional | container recipe, if your method needs its own image |

## 2. The rule contract (`run_<method>.smk`)

Minimal R example, modeled on `nnls`:

```python
import os, yaml
with open("config.yaml") as f:
    params = yaml.safe_load(f)

def get_basename(p):
    return os.path.splitext(os.path.basename(p))[0]

sc_input    = config["sc_input"]
sp_input    = config["sp_input"]
output_dir  = config["output"]
runID_props = params["runID_props"]
method      = "<method>"
output      = f"{output_dir}/proportions_{method}_{get_basename(sp_input)}{runID_props}.tsv"
annot       = config.get("annot", params["annot"])
map_genes   = config.get("map_genes", "false")

rule run_<method>:
    input:
        sc_input=sc_input,
        sp_input=sp_input
    output:
        output
    singularity:
        "docker://ghcr.io/cbib/sp_<method>:latest"   # any docker:// image
    threads: 8
    shell:
        """
        Rscript subworkflows/deconvolution/<method>/script.R \
            --sc_input {input.sc_input} --sp_input {input.sp_input} \
            --annot {annot} --output {output} --map_genes {map_genes}
        """
```

**Contract — respect these or the pipeline won't find your results:**

- **Output path is fixed**: `<output>/proportions_<method>_<sample>{runID_props}.tsv`. `main.smk`
  builds exactly this path for every requested method — the rule's `output:` must match it.
- **Output format**: TSV, **spots in rows, cell types in columns**, each row summing to 1. Copy the
  result-writing block from an existing `script.R` (it omits rownames quoting, strips
  non-alphanumeric characters from cell-type names, and shell-sorts the columns) so every method's
  columns line up in the visualization.
- **`method.env`**: `USE_GPU=true` → GPU partition (Singularity runs with `--nv`); `false` → CPU.

## 3. Your script

The rule passes your script:

- `--sc_input` — single-cell reference `.rds` (a Seurat object with a cell-type annotation column).
- `--sp_input` — spatial `.rds` (a Seurat object, **or** a synthspot named list with the expression
  matrix in `counts`).
- `--annot` — name of the cell-type column in the reference metadata.
- `--output` — where to write the TSV.
- `--map_genes` — `true` maps gene symbols ↔ Ensembl when the reference and spatial data use
  different gene IDs. The cbib R images bundle `org.Hs.eg.db` so this works offline; a Python
  method can use `pybiomart` (as `cell2location` does).

## 4. Python methods with a build / fit split

Bayesian methods (e.g. `cell2location`) separate model **building** (single-cell only) from
**fitting** (spatial), so the model is built once and reused across spatial samples. Write two
rules in `run_<method>.smk` chained by their outputs — see
`subworkflows/deconvolution/cell2location/run_cell2location.smk` as the reference implementation
(note it invokes the container's env Python explicitly, e.g.
`/opt/conda/envs/<env>/bin/python`, not bare `python3`).

## 5. Test it

```bash
# via the launcher (routes CPU/GPU from method.env), running locally:
./deconvolista --local --methods <method> \
  --sc unit-test/test_sc_data.rds --sp unit-test/test_sp_data.rds --annot subclass

# or plain Snakemake:
snakemake -s main.smk -c4 --use-singularity \
  --config mode=run_dataset methods=<method> \
  sc_input=unit-test/test_sc_data.rds sp_input=unit-test/test_sp_data.rds \
  output=res annot=subclass
```

Then `--methods <method>` is immediately usable everywhere — no other file to touch.

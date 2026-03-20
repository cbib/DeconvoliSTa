# Rule: run_dataset

## Purpose  
Performs deconvolution on a spatial transcriptomics dataset. For every method listed in `methods`, a proportion table is produced.  
If the global option `skip_metrics` is **false**, evaluation rules are automatically added.

### Required configuration (config.yaml or `--config`)
| Key          | Type    | Example                              |
|--------------|---------|--------------------------------------|
| mode         | string  | `"run_dataset"`                      |
| sc_input     | file    | `"reference_sc.rds"`                 |
| sp_input     | file    | `"sample_sp.rds"`                    |
| methods      | csv list | `"cell2location,rctd,nnls"`          |

### Optional configuration
| Key            | Default | Notes                                            |
|----------------|---------|--------------------------------------------------|
| skip_metrics   | `"false"` | Set to `"true"` to skip evaluation metrics.    |
| output         | `.`     | Directory where `.tsv` files are written.       |
| runID_props    | ``      | Optional suffix for filenames (e.g., “_v1”).    |

### Example
```bash
snakemake -s run_dataset.smk \
  --config mode="run_dataset" methods=cell2location,rctd,nnls \
  sc_input="ref_sc.rds" sp_input="sample_sp.rds" output="results"
```

---


# Rule: generate_data

## Purpose  
Creates synthetic spatial datasets from a single‑cell reference using SynthSpot.  

### Required configuration
| Key          | Type    | Example                     |
|--------------|---------|-----------------------------|
| mode         | string  | `"generate_data"`           |
| sc_input     | file    | `"gold_standard_1.rds"`     |
| dataset_type | csv list| `"aud,add"`                 |
| reps         | int     | `5`                         |
| rootdir      | dir     | `"synthetic_data_sm"`       |

> **Dataset type mapping** – abbreviations map to full names as defined in the Snakefile (`synthspot_types_map`).  

### Optional configuration
| Key   | Default | Notes                            |
|-------|---------|----------------------------------|
| output| `.`     | Directory where synthetic files are stored. |

### Example
```bash
snakemake -s generate_data.smk \
  --config mode="generate_data" sc_input="gold_standard_1.rds" dataset_type=aud,add reps=3 rootdir="synthetic_data_sm"
```

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

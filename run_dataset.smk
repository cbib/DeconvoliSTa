import os, sys, yaml

# --- shared helper functions -------------------------------------------------
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

def get_config_var(config, var_name, default=None):
    if var_name not in config:
        if default is None:
            raise ValueError(f"Error: '{var_name}' is missing in the configuration file.")
        print(f"Warning: '{var_name}' is missing. Using default: {default}")
        return default
    return config[var_name]

# --- load configuration ------------------------------------------------------
with open("config.yaml", "r") as f:
    cfg = yaml.safe_load(f)

mode      = get_config_var(cfg, "mode")
sc_input  = get_config_var(cfg, "sc_input")
sp_input  = get_config_var(cfg, "sp_input")
methods   = get_config_var(cfg, "methods").split(',')
supported_methods = {"cell2location","rctd","spatialdwls","nnls","ddls","dirichlet","seurat"}
unsupported = set(methods) - supported_methods
if unsupported:
    raise ValueError(f"Unsupported methods: {', '.join(unsupported)}")

output_dir   = get_config_var(cfg, "output")
skip_metrics = get_config_var(cfg, "skip_metrics", "false")
runID_props  = get_config_var(cfg.get("params", {}), "runID_props")

# --- file lists ---------------------------------------------------------------
output_suffix = get_basename(sp_input)
proportions_files = [
    f"{output_dir}/proportions_{m}_{output_suffix}{runID_props}.tsv" for m in methods
]

if skip_metrics == "false":
    metrics_files = [
        f"{output_dir}/metrics/metrics_{m}_{output_suffix}{runID_props}.tsv"
        for m in methods
    ]

# --- rules -------------------------------------------------------------------
rule all:
    input:
        *proportions_files,
        *metrics_files if skip_metrics == "false" else []

include: "subworkflows/deconvolution/run_methods.smk"

if skip_metrics == "false":
    include: "subworkflows/evaluation/evaluate_methods.smk"
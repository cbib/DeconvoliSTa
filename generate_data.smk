import os, sys, yaml

def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

def get_config_var(config, var_name, default=None):
    if var_name not in config:
        if default is None:
            raise ValueError(f"Error: '{var_name}' missing.")
        print(f"Warning: '{var_name}' missing. Using {default}")
        return default
    return config[var_name]

Need to add all values into config.yaml? Actually rewrite code to read config with full mapping. Replace open line with reading entire file and then load yaml. Provide code snippet.```smk
import yaml

# Load configuration from YAML file
with open("config.yaml", "r") as f:
    cfg = yaml.safe_load(f)

# Mapping of dataset type abbreviations to full names
synthspot_types_map = {
    'aud': "artificial_uniform_distinct",
    'add': "artificial_diverse_distinct",
    'auo': "artificial_uniform_overlap",
    'ado': "artificial_diverse_overlap",
    'adcd': "artificial_dominant_celltype_diverse",
    'apdcd': "artificial_partially_dominant_celltype_diverse",
    'adrcd': "artificial_dominant_rare_celltype_diverse",
    'arrcd': "artificial_regional_rare_celltype_diverse",
    'prior': "prior_from_data",
    'real': "real",
    'rm': "real_missing_celltypes_visium",
    'arm': "artificial_missing_celltypes_visium",
    'addm': "artificial_diverse_distinct_missing_celltype_sc",
    'adom': "artificial_diverse_overlap_missing_celltype_sc"
}
    cfg = yaml.safe_load(f)

 

mode       = get_config_var(cfg, "mode")
sc_input   = get_config_var(cfg, "sc_input")
dataset_type_str = get_config_var(cfg, "dataset_type")
dataset_types = [synthspot_types_map[t] for t in dataset_type_str.split(',')]
reps        = int(get_config_var(cfg, "reps"))
rootdir     = get_config_var(cfg, "rootdir")

generated_files = [
    f"synthetic_data_sm/{os.path.basename(sc_input).split('.')[0]}_{dt}_rep{r}.rds"
    for dt in dataset_types
    for r in range(1, reps+1)
]

include: "subworkflows/data_generation/generate_data.smk"

rule all:
    input: generated_files
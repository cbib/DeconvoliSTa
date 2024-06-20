import os
import sys
import yaml

# Lire le fichier de configuration YAML
with open("my_config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)
# Fonction pour obtenir le nom de base du fichier sans extension
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

# Charger les paramètres de configuration nécessaires
sc_input = config["sc_input"]
mode = config["mode"]


if mode  == "run_dataset":
    sp_input = config["sp_input"]
    methods =  config["methods"].split(',')   

    output_suffix = get_basename(sp_input)
    runID_props = params["runID_props"]
    include: "subworkflows_sm/deconvolution/run_methods.smk"
    include: "subworkflows_sm/evaluation/evaluate_methods.smk"
    output_dir = config["output"]

    output= [f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv" for method in methods]
    metrics_files = [f"{output_dir}/metrics/metrics_{method}_{output_suffix}{runID_props}.tsv" for method in methods]

    rule main:
        input:
            output_files,
            metrics_files
elif mode == "generate_data":
    generated_files = config["sc_input"]
    dataset_types = config["dataset_type"].split(',')
    reps = config["reps"]
    rootdir = config["rootdir"]

    generated_files = [f"synthetic_data_sm/{os.path.basename(sc_input).split('.')[0]}_{dataset_type}_rep{rep}.rds" for dataset_type in dataset_types for rep in range(1, int(reps) + 1)]
    include: "subworkflows_sm/data_generation/generate_data.smk"
    rule gen_files:
        input:
            generated_files
else:
    print("Enter a valid execution mode --mode\n")

import os
import sys
import yaml

# Lire le fichier de configuration YAML
with open("subworkflows_sm/subworkflows/deconvolution/cell2location/my_config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)
# Fonction pour obtenir le nom de base du fichier sans extension
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

# Charger les paramètres de configuration nécessaires
methods =  config["methods"].split(',')   
sc_input = config["sc_input"]
sp_input = config["sp_input"]

output_suffix = get_basename(sp_input)
runID_props = params["runID_props"]


include: "subworkflows_sm/subworkflows/deconvolution/run_methods.smk"
output_dir = config["output"]


output= [f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv" for method in methods]

rule main:
    input:
        output

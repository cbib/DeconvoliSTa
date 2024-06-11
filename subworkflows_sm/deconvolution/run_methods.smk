
import os
import sys
import yaml

# Lire le fichier de configuration YAML
with open("subworkflows_sm/deconvolution/cell2location/my_config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

# Fonction pour obtenir le nom de base du fichier sans extension
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

# Charger les paramètres de configuration nécessaires
methods = config["methods"].split(',')   
sc_input = config["sc_input"]
sp_input = config["sp_input"]
output_dir = config["output"]

output_suffix = get_basename(sp_input)
runID_props = params["runID_props"]

# Inclure dynamiquement les pipelines appropriés en fonction des méthodes
for method in methods:
    method = method.strip()  
    include: f"{method}/run_{method}.smk"

# Générer une liste des fichiers de sortie pour chaque méthode
# output_files = [f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.preformat" for method in methods]
output_files = [f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv" for method in methods]

for f in output_files:
    rule_name = f
    rule rule_name:
        input:
            f

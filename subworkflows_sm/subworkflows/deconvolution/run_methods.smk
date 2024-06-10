import os
import sys
import yaml

# Lire le fichier de configuration YAML
with open("cell2location/my_config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

# Charger les paramètres de configuration nécessaires
methods =  config["methods"] else "cell2location"  # Par défaut, utiliser "cell2location" si non spécifié
sc_input = config["sc_input"]
sp_input = config["sp_input"]
output_suffix = get_basename(sp_input)
runID_props = params["runID_props"]

# Définir la règle principale qui inclut le pipeline approprié en fonction de la méthode
if methods == "cell2location":
    include: "cell2location/run_cell2location.smk"
elif methods == "RCTD":
    include: "rctd/run_rctd.smk"
else:
    raise ValueError(f"Unsupported methods: {methods}")

# Fonction pour obtenir le nom de base du fichier sans extension
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

output= f"proportions_{methods}_{output_suffix}{runID_props}.preformat"

# Règle all générique qui pourrait être adaptée en fonction de la méthode
rule all:
    input:
        out = output

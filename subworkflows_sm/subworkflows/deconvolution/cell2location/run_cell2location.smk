# Snakefile

import os
import sys
# sys.path.append('.')
from functions import build_cell2location_model, fit_cell2location_model


# # Préparez les chemins d'entrée/sortie
sc_input = config["sc_input"]
sp_input = config["sp_input"]
output = config["output"]
import yaml

# Lire le fichier de configuration YAML
with open("my_config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

# rule all:
#     input:
#         "sc.h5ad",
#         # "proportions_cell2location_{output_suffix}{runID_props}.preformat"

# Fonction pour obtenir le nom de base du fichier sans extension
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

rule convertBetweenRDSandH5AD:
    input:
        rds_file=sc_input
    output:
        h5ad_file= f"{get_basename(config["sc_input"])}.h5ad",
    singularity:
        "docker://csangara/seuratdisk:latest"
    script:
        """
        Rscript ./subworkflows/deconvolution/convertBetweenRDSandH5AD.R --input_path {input.rds_file}
        """

# rule build_cell2location:
#     input:
#         sc_input
#     output:
#         "sc.h5ad"
#     singularity:
#         "docker://csangara/sp_cell2location:latest"
#     shell:
#         """
#         python3 functions.py {sc_input} 
#         """


# rule fit_cell2location:
#     input:
#         sp_input,
#         model="sc.h5ad"
#     output:
#         "proportions_cell2location_{output_suffix}{runID_props}.preformat"
#     singularity:
#         "docker://csangara/sp_cell2location:latest"
#     shell:
#         """
#         fit_cell2location_model {sp_input}, {model}, {config}
#         """


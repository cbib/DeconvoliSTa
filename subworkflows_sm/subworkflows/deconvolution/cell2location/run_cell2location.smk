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

# Fonction pour convertir un fichier RDS en H5AD
def convert_rds_to_h5ad(rds_file, h5ad_file):
    import rpy2.robjects as robjects
    from rpy2.robjects import r
    r['load'](rds_file)
    adata = robjects.r('data')  # Assuming the loaded object is named 'data'
    adata.write_h5ad(h5ad_file)

# Vérifiez si le fichier sc_input est en format RDS et convertissez-le en H5AD si nécessaire
if sc_input.endswith('.rds'):
    h5ad_file = sc_input.replace('.rds', '.h5ad')
    convert_rds_to_h5ad(sc_input, h5ad_file)
    sc_input = h5ad_file

rule build_cell2location:
    input:
        sc_input
    output:
        "sc.h5ad"
    singularity:
        "docker://csangara/sp_cell2location:latest"
    shell:
        """
        python3 functions.py {sc_input} 
        """


rule fit_cell2location:
    input:
        sp_input,
        model="sc.h5ad"
    output:
        "proportions_cell2location_{output_suffix}{runID_props}.preformat"
    singularity:
        "docker://csangara/sp_cell2location:latest"
    shell:
        """
        fit_cell2location_model {sp_input}, {model}, {config}
        """

# rule format_c2l:
#     input:
#         "proportions_cell2location_{output_suffix}{runID_props}.preformat"
#     output:
#         "formatted_proportions_cell2location_{output_suffix}{runID_props}.tsv"
#     shell:
#         format_tsv(input, output, config)



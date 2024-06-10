# Snakefile

import os
import sys
# sys.path.append('.')
from run_build import build_cell2location_model 
from run_fit import fit_cell2location_model

# # Préparez les chemins d'entrée/sortie
sc_input = config["sc_input"]
sp_input = config["sp_input"]
output = config["output"]
import yaml

# Lire le fichier de configuration YAML
with open("my_config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)


# Fonction pour obtenir le nom de base du fichier sans extension
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

output_suffix = get_basename(sp_input)
runID_props = params["runID_props"]

rule all:
    input:
        "proportions_cell2location_{}{}.preformat".format(output_suffix, runID_props)
    run:
        print(f"output_suffix: {output_suffix}")
        print(f"runID_props: {runID_props}")

rule convertBetweenRDSandH5AD:
    input:
        sc_rds_file=sc_input,
        sp_rds_file=sp_input
    output:
        sc_h5ad_file=temp(f"{get_basename(sc_input)}.h5ad"),
        sp_h5ad_file=temp(f"{get_basename(sp_input)}.h5ad")
    singularity:
        "docker://csangara/seuratdisk:latest"
    shell:
        r"""
        Rscript ./convertBetweenRDSandH5AD.R --input_path {input.sc_rds_file} 
        """
        r"""
        Rscript ./convertBetweenRDSandH5AD.R --input_path {input.sp_rds_file} 
        """
rule build_cell2location:
    input:
        rules.convertBetweenRDSandH5AD.output.sc_h5ad_file
    output:
        "sc.h5ad"
    singularity:
        "docker://csangara/sp_cell2location:latest"
    shell:
        """
        python3 run_build.py {input[0]}
        """


rule fit_cell2location:
    input:
        rules.convertBetweenRDSandH5AD.output.sp_h5ad_file,
        model="sc.h5ad"
    output:
        "proportions_cell2location_{}{}.preformat".format(output_suffix, runID_props)
    singularity:
        "docker://csangara/sp_cell2location:latest"
    shell:
        """
        python3 run_fit.py {input[0]} {input[1]}
        """


# Snakefile

import os
import yaml
import time
# Lire le fichier de configuration YAML
with open("config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

# Fonction pour obtenir le nom de base du fichier sans extension
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]


# Helper function to check if a config variable is present
def get_config_var(config, var_name, default=None):
    if var_name not in config:
        if default is None:
            raise ValueError(f"Error: '{var_name}' is missing in the configuration file.")
        else:
            print(f"Warning: '{var_name}' is missing in the configuration file. Using default: {default}")
            return default
    return config[var_name]

# Charger les paramètres de configuration nécessaires
sc_input = get_config_var(config, "sc_input")
sp_input = get_config_var(config, "sp_input")
output_dir = get_config_var(config, "output")
output_suffix = get_basename(sp_input)
runID_props = get_config_var(params, "runID_props")
method = "spatialdwls"
output = f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv"
deconv_args = get_config_var(params, "deconv_args")

# Définir le chemin absolu du script R
script_dir = os.path.dirname(os.path.abspath(__file__))
rctd_script = "subworkflows/deconvolution/spatialdwls/script_nf.R"

annot = get_config_var(config, "annot", get_config_var(params, "annot"))
map_genes = get_config_var(config, "map_genes", "false")

rule run_spatialdwls:
    input:
        sc_input=sc_input,
        sp_input=sp_input
    output:
        output
    singularity:
        "docker://csangara/sp_spatialdwls:latest" #"docker://csangara/sp_rctd:latest"
    threads:
        12
    shell:
        """
        Rscript subworkflows/deconvolution/spatialdwls/script_nf.R \
            --sc_input {sc_input} --sp_input {sp_input} \
            --annot {annot} --output {output}  --map_genes {map_genes}
        """


# Snakefile
import os
import yaml

with open("my_config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

# # Fonction pour obtenir le nom de base du fichier sans extension
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

sc_input = config["sc_input"]
sp_input = config["sp_input"]
output_dir = config["output"]
output_suffix = get_basename(sp_input)
runID_props = params["runID_props"]
method = "dirichlet"
formatted_output = f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv"
use_gpu = config["use_gpu"]

# Définir le chemin absolu du script R
script_dir = os.path.dirname(os.path.abspath(__file__))
convert_script = "subworkflows/deconvolution/convertBetweenRDSandH5AD.R"
dirichlet_script = "subworkflows/deconvolution/dirichlet/gen_dirichlet.R"

rule run_dirichlet:
    input:
        sp_input = sp_input
    output:
        formatted_output
    singularity:
        "docker://csangara/seuratdisk:latest"
    shell:
        """
        Rscript {dirichlet_script} --sp_input {input.sp_input} --output {formatted_output}
        """

# Snakefile

import os
import yaml

# Lire le fichier de configuration YAML
with open("subworkflows_sm/deconvolution/cell2location/my_config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

# Fonction pour obtenir le nom de base du fichier sans extension
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

sc_input = config["sc_input"]
sp_input = config["sp_input"]
output_dir = config["output"]
output_suffix = get_basename(sp_input)
runID_props = params["runID_props"]
method = "nnls"
output = f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv"
deconv_args = params['deconv_args']
# Définir le chemin absolu du script R
script_dir = os.path.dirname(os.path.abspath(__file__))
nnls_script = "subworkflows_sm/deconvolution/nnls/script_nf.R"
annot = config["annot"] if "annot" in config.keys() else params['annot']
rule nnls:
    input:
        output

rule run_nnls:
    input:
        sc_input=sc_input,
        sp_input=sp_input
    output:
        output
    singularity:
        "docker://csangara/sp_nnls:latest"
    threads:
        8 
    shell:
        """
        Rscript {nnls_script} \
             --sc_input {input.sc_input} --sp_input {input.sp_input} \
            --annot {annot} --output {output}
        """

# Snakefile

import os
import yaml
import time
# Read the YAML configuration file
with open("config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

# Function to get the file basename (without extension)
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

sc_input = config["sc_input"]
sp_input = config["sp_input"]
output_dir = config["output"]
output_suffix = get_basename(sp_input)
runID_props = params["runID_props"]
method = "ddls"
output = f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv"
deconv_args = params['deconv_args']
# Absolute path to the R script
ddls_script = "subworkflows/deconvolution/ddls/script.R"
annot = config["annot"] if "annot" in config.keys() else params['annot']
use_gpu = config["use_gpu"] if "use_gpu" in config.keys() else "false"

rule run_ddls:
    input:
        sc_input=sc_input,
        sp_input=sp_input
    output:
        output
    singularity:
        "docker://abderahim02/sp_ddls:latest"
    threads:
        12
    shell: # --output {output} --epochs 10000 --batch_size 10
        """
        Rscript {ddls_script} \
            --sc_input {input.sc_input} --sp_input {input.sp_input} \
            --output {output} --epochs 20000 --batch_size 15 --use_gpu {use_gpu}
        """
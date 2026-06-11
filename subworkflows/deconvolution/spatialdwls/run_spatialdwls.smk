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


# Helper function to check if a config variable is present
def get_config_var(config, var_name, default=None):
    if var_name not in config:
        if default is None:
            raise ValueError(f"Error: '{var_name}' is missing in the configuration file.")
        else:
            print(f"Warning: '{var_name}' is missing in the configuration file. Using default: {default}")
            return default
    return config[var_name]

# Load the required configuration parameters
sc_input = get_config_var(config, "sc_input")
sp_input = get_config_var(config, "sp_input")
output_dir = get_config_var(config, "output")
output_suffix = get_basename(sp_input)
runID_props = get_config_var(params, "runID_props")
method = "spatialdwls"
output = f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv"
deconv_args = get_config_var(params, "deconv_args")

# Absolute path to the R script
script_dir = os.path.dirname(os.path.abspath(__file__))
spatialdwls_script = "subworkflows/deconvolution/spatialdwls/script.R"

annot = get_config_var(config, "annot", get_config_var(params, "annot"))
map_genes = get_config_var(config, "map_genes", "false")

# Use a locally built image if it exists in SIF_DIR, otherwise pull the public one.
# NB: gene mapping (map_genes=true) needs the org.Hs.eg.db-enabled local build (spatialdwls.def);
# the public image works for data whose genes already match the reference.
_local_sif = f"{config.get('sif_dir', 'sif')}/sp_spatialdwls_cbib.sif"
spatialdwls_image = _local_sif if os.path.exists(_local_sif) else "docker://csangara/sp_spatialdwls:latest"

rule run_spatialdwls:
    input:
        sc_input=sc_input,
        sp_input=sp_input
    output:
        output
    singularity:
        spatialdwls_image
    threads:
        12
    shell:
        """
        Rscript {spatialdwls_script} \
            --sc_input {input.sc_input} --sp_input {input.sp_input} \
            --annot {annot} --output {output}  --map_genes {map_genes}
        """


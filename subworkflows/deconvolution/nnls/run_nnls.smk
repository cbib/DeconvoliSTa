# Snakefile

import os
import yaml

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
method = "nnls"
output = f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv"
deconv_args = params['deconv_args']
# Absolute path to the R script
script_dir = os.path.dirname(os.path.abspath(__file__))
nnls_script = "subworkflows/deconvolution/nnls/script.R"
annot = config["annot"] if "annot" in config.keys() else params['annot']
map_genes = config.get("map_genes", "false")

# cbib image on GHCR (csangara + org.Hs.eg.db, so map_genes=true works out of the box).
# NB: ghcr-only for now to test the registry; the local-sif override will be re-added later.
nnls_image = "docker://ghcr.io/cbib/sp_nnls:latest"

rule run_nnls:
    input:
        sc_input=sc_input,
        sp_input=sp_input
    output:
        output
    singularity:
        nnls_image
    threads:
        8 
    shell:
        """
        Rscript {nnls_script} \
             --sc_input {input.sc_input} --sp_input {input.sp_input} \
            --annot {annot} --output {output} --map_genes {map_genes}
        """

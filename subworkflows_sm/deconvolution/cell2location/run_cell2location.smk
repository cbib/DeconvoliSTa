import os
import yaml
import time

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
method = "cell2location"
formatted_output = f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv"
use_gpu = config["use_gpu"]

# Définir le chemin absolu du script R
script_dir = os.path.dirname(os.path.abspath(__file__))
convert_script = "subworkflows_sm/deconvolution/convertBetweenRDSandH5AD.R"

rule all:
    input:
        formatted_output

rule convertBetweenRDSandH5AD:
    input:
        sc_rds_file=sc_input,
        sp_rds_file=sp_input
    output:
        sc_h5ad_file=temp(f"{get_basename(sc_input)}.h5ad"),
        sp_h5ad_file=temp(f"{get_basename(sp_input)}.h5ad")
    singularity:
        "docker://csangara/seuratdisk:latest"
    threads:
        8
    shell:
        """
        start_time=$(date +%s)
        Rscript {convert_script} --input_path {input.sc_rds_file}
        Rscript {convert_script} --input_path {input.sp_rds_file}
        end_time=$(date +%s)
        elapsed_time=$((end_time - start_time))
        echo "convertBetweenRDSandH5AD took $elapsed_time seconds"
        """

rule build_cell2location:
    input:
        rules.convertBetweenRDSandH5AD.output.sc_h5ad_file
    output:
        temp(f"{output_dir}/sc.h5ad")
    singularity:
        "docker://csangara/sp_cell2location:latest"
    threads:
        8
    shell:
        """
        start_time=$(date +%s)
        python3 subworkflows_sm/deconvolution/cell2location/run_build.py {input[0]} {output_dir} {use_gpu}
        end_time=$(date +%s)
        elapsed_time=$((end_time - start_time))
        echo "build_cell2location took $elapsed_time seconds"
        """

rule fit_cell2location:
    input:
        rules.convertBetweenRDSandH5AD.output.sp_h5ad_file,
        model=f"{output_dir}/sc.h5ad"
    output:
        temp(f"{output_dir}/proportions_cell2location_{output_suffix}{runID_props}.preformat")
    singularity:
        "docker://csangara/sp_cell2location:latest"
    threads:
        8
    shell:
        """
        start_time=$(date +%s)
        python3 subworkflows_sm/deconvolution/cell2location/run_fit.py {input[0]} {input[1]} {output_dir} {use_gpu}
        end_time=$(date +%s)
        elapsed_time=$((end_time - start_time))
        echo "fit_cell2location took $elapsed_time seconds"
        """

rule format_tsv_file:
    input:
        tsv_file=rules.fit_cell2location.output
    output:
        formatted_output
    singularity:
        "docker://rocker/tidyverse:latest"
    threads:
        8
    shell:
        """
        start_time=$(date +%s)
        Rscript -e "
        deconv_matrix <- read.table('{input.tsv_file}', sep='\t', header=TRUE, row.names=1);
        colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), '[/. ]', '');
        deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method='shell')];
        write.table(deconv_matrix, file='{output}', sep='\t', quote=FALSE, row.names=FALSE);
        "
        end_time=$(date +%s)
        elapsed_time=$((end_time - start_time))
        echo "format_tsv_file took $elapsed_time seconds "
        """

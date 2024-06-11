import os
import yaml

# Charger le fichier de configuration
with open("subworkflows_sm/subworkflows/evaluation/my_config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

# Préparez les chemins d'entrée/sortie
sc_input = config["sc_input"]
sp_input = config["sp_input"]
output_dir = config["output"]
methods =  config["methods"].split(',')   
runID_props = params["runID_props"]
runID_metrics = params["runID_metrics"]
# p = params["params"]
remap_annot = params["params"]["remap_annot"]
# Fonction pour obtenir le nom de base du fichier sans extension
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

output_suffix = get_basename(sp_input)
print("methods evaluate = ", methods)
# Générer les fichiers de métriques pour chaque méthode
metrics_files = [f"{output_dir}/metrics/metrics_{method}_{output_suffix}{runID_props}.tsv" for method in methods]
output_files= [f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv" for method in methods]

rule all_metrics:
    input:
        output_files, 
        metrics_files

# Générer dynamiquement les règles pour chaque méthode
for method in methods:
    rule:
        input:
            props_file=f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv",
            sp_input=sp_input
        output:
            metrics_file=f"{output_dir}/metrics/metrics_{method}_{output_suffix}{runID_props}.tsv"
        singularity:
            "docker://csangara/sp_eval:latest"
        shell:
            """
            Rscript subworkflows_sm/subworkflows/evaluation/metrics.R \
            --props_file {input.props_file} --sp_input {input.sp_input} \
            --output {output.metrics_file}
            """

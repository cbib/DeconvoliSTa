# # Snakefile
# import numpy as np
# import pandas as pd

# def generate_dirichlet_dataframe(alpha, rows, columns):
#     """
#     Genere un DataFrame de valeurs suivant une distribution de Dirichlet.

#     :param alpha: Parametre de la distribution de Dirichlet. Peut etre un seul float ou une liste de floats.
#     :param rows: Nombre de lignes du DataFrame.
#     :param columns: Liste des noms de colonnes du DataFrame.
#     :return: Un DataFrame Pandas avec les valeurs generees selon la distribution de Dirichlet.
#     """
#     # Vérifier que la longueur de alpha correspond au nombre de colonnes
#     if isinstance(alpha, (list, np.ndarray)):
#         assert len(alpha) == len(columns), "La longueur de alpha doit correspondre au nombre de colonnes."
#     dirichlet_values = np.random.dirichlet(alpha, size=rows)
#     df = pd.DataFrame(dirichlet_values, columns=columns)
#     return df

# def save_dataframe_to_tsv(df, file_path):
#     """
#     Sauvegarde un DataFrame Pandas dans un fichier TSV.
#     """
#     df.to_csv(file_path, sep='\t', index=False)


import os
# import yaml

# Lire le fichier de configuration YAML
# with open("subworkflows_sm/deconvolution/cell2location/my_config.yaml", "r") as config_file:
#     params = yaml.safe_load(config_file)

# # Fonction pour obtenir le nom de base du fichier sans extension
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

# sc_input = config["sc_input"]
sp_input = config["sp_input"]
# output_dir = config["output"]
# output_suffix = get_basename(sp_input)
# runID_props = params["runID_props"]
# method = "dirichlet"
# formatted_output = f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv"
# use_gpu = config["use_gpu"]

# Définir le chemin absolu du script R
script_dir = os.path.dirname(os.path.abspath(__file__))
convert_script = "subworkflows_sm/deconvolution/convertBetweenRDSandH5AD.R"

# rule all:
#     input:
#         formatted_output

rule convertBetweenRDSandH5AD:
    input:
        sp_rds_file=sp_input
    output:
        sp_h5ad_file=f"{get_basename(sp_input)}.h5ad"
    singularity:
        "docker://csangara/seuratdisk:latest"
    shell:
        """
        Rscript {convert_script} --input_path {input.sp_rds_file}
        """

# rule run_dirichlet:
#     input:
#         rules.convertBetweenRDSandH5AD.output
#     output:
#         formatted_output
#     run:
#         alpha = [0.5, 0.5, 0.5]  
#         rows = 5  
#         columns = ['A', 'B', 'C'] 
#         df = generate_dirichlet_dataframe(alpha, rows, columns)
#         save_dataframe_to_tsv(df, file_path)



# functions.py
import subprocess
import configparser
import yaml

# Fonction pour obtenir le nom de base du fichier sans extension
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

# # Lire le fichier de configuration YAML
with open("subworkflows_sm/deconvolution/cell2location/my_config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

def fit_cell2location_model(sp_input, model, output_dir, use_gpu, map_genes):
    """
    Fit cell2location model.
    
    Args:
        sp_input (str): Path to spatial input file.
        model (str): Path to the model file.
        params (dict): Parameters for fitting the model.
    """
    output_suffix = get_basename(sp_input)
    output = f"proportions_cell2location_{output_suffix}{params['runID_props']}.preformat"
    epochs = f"-e {params['epoch_fit']}" if params['epoch_fit'] != "default" else ""
    args = params.get('deconv_args', {}).get('cell2location', {}).get('fit', "")
    # p_parameter = f"-p 5"

    cuda_device = params["cuda_device"] if use_gpu == "true" else "cpu"
    # output_dir = params.get('output_dir', '.')
    
    run_dev = 'GPU' if use_gpu == "true" else 'CPU'
    print(f"Fitting cell2location model from file {model} with {run_dev}...")
    print(f"Arguments: {args}")
    print(f"{sp_input}")
    # print(f"sp_input.split(',')[0]" {sp_input.split(",")[0]})
    command = [ #pip install pybiomart 
        "bash", "-c", f"source activate cell2loc_env && python subworkflows_sm/deconvolution/cell2location/fit_model.py {sp_input.split(',')[0]} {model} {cuda_device} {epochs} {args} -o {output_dir}  -m {map_genes} && mv {output_dir}/proportions.tsv {output_dir}/{output}"
    ]
    print(command)
    subprocess.run(command, check=True)

def format_tsv(input_file, output_file, params):
    # Assuming you have a script or function for formatting
    command = f"python format_script.py {input_file} -o {output_file}"
    subprocess.run(command, shell=True, check=True)

if __name__ == "__main__":
    import os
    import sys  
    # Récupérer tous les arguments de la ligne de commande
    args = sys.argv
    sp_input  = args[1]
    model = args[2]
    output_dir = args[3]
    use_gpu = args[4]
    map_genes = args[5]
    fit_cell2location_model(sp_input, model, output_dir, use_gpu, map_genes)
  
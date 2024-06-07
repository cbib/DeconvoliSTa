# functions.py
import subprocess
import configparser
import yaml

# # Lire le fichier de configuration YAML
with open("my_config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

def fit_cell2location_model(sp_input, model):
    """
    Fit cell2location model.
    
    Args:
        sp_input (str): Path to spatial input file.
        model (str): Path to the model file.
        params (dict): Parameters for fitting the model.
    """
    output_suffix = sp_input.split('/')[-1]
    output = f"proportions_cell2location_{output_suffix}{params['runID_props']}.preformat"
    epochs = f"-e {params['epoch_fit']}" if params['epoch_fit'] != "default" else ""
    args = params.get('deconv_args', {}).get('cell2location', {}).get('fit', "")
    cuda_device = "cpu"
    output_dir = params.get('output_dir', '.')
    
    print(f"Fitting cell2location model from file {model} with {'GPU' if params['gpu'] else 'CPU'}...")
    print(f"Arguments: {args}")
    
    command = [
        "bash", "-c", f"source activate cell2loc_env && python fit_model.py {sp_input} {model} {cuda_device} {epochs} {args} -o {output_dir} -p 5 && mv {output_dir}/proportions.tsv {output}"
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
    sc_input  = args[1]
    model = args[2]
    fit_cell2location_model(sc_input, model)
  
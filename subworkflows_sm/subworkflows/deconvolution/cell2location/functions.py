# functions.py

import subprocess

import configparser

# Lire le fichier de configuration Nextflow-like
# conf = configparser.ConfigParser()
# conf.read("../my_config.config")
# # params = conf["params"]

import yaml

# Lire le fichier de configuration YAML
with open("my_config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

def build_cell2location_model(sc_input):
    """
    Build cell2location model.
    
    Args:
        sc_input (str): Path to single-cell input file.
        params (dict): Parameters for building the model.
    """
    tag_suffix = sc_input.split('/')[-1]
    sample_id_arg = f"-s {params['sampleID']}" if params['sampleID'] != "none" else ""
    epochs = f"-e {params['epoch_build']}" if params['epoch_build'] != "default" else ""
    args = params.get('deconv_args', {}).get('cell2location', {}).get('build', "")
    cuda_device = "cpu"
    output_dir = params.get('output_dir', '.')
    
    print(f"Building cell2location model with {'GPU' if params['gpu'] else 'CPU'}...")
    
    command = [
        "bash", "-c", f"""
        source activate cell2loc_env &&
        python build_model.py -a {params['annot']} {sample_id_arg} {epochs} {args} -o {output_dir}
            {sc_input} {cuda_device} 
        """
    ]
    print(command)
    subprocess.run(command, check=True)

def fit_cell2location_model(sp_input, model, params):
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
        "bash", "-c", f"""
        source activate cell2loc_env &&
        python fit_model.py \\
            {sp_input} {model} {cuda_device} {epochs} {args} -o {output_dir} &&
        mv {output_dir}/proportions.tsv {output}
        """
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
    print(sc_input)

    build_cell2location_model(sc_input)
  
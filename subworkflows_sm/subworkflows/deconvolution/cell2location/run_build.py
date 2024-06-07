# functions.py
import subprocess
import configparser
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
    import os
    command = [
        "bash", "-c", f"source activate cell2loc_env && python build_model.py {sc_input} {cuda_device} -a {params['annot']} {sample_id_arg} {epochs} {args} -o {output_dir} -p 5"
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
  
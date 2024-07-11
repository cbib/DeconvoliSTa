# functions.py
import subprocess
import configparser
import yaml


# Lire le fichier de configuration YAML
with open("subworkflows_sm/deconvolution/cell2location/my_config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

def build_cell2location_model(sc_input, output_dir, use_gpu, annot):
    """
    Build cell2location model.
    
    Args:
        sc_input (str): Path to single-cell input file.
    """
    tag_suffix = sc_input.split('/')[-1]
    sample_id_arg = f"-s {params['sampleID']}" if params['sampleID'] != "none" else ""
    epochs = f"-e {params['epoch_build']}" if params['epoch_build'] != "default" else ""
    args = params.get('deconv_args', {}).get('cell2location', {}).get('build', "")
    cuda_device = params["cuda_device"] if use_gpu == "true" else "cpu"
    # p_parameter = f"-p 5"
    run_dev = 'GPU' if use_gpu == "true" else 'CPU'
    print(f"Building cell2location model with {run_dev}...")
    import os
    command = [
        "bash", "-c", f"source activate cell2loc_env && python subworkflows_sm/deconvolution/cell2location/build_model.py {sc_input} {cuda_device} -a {annot} {sample_id_arg} {epochs} {args} -o {output_dir} "
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
    output_dir = args[2]
    use_gpu = args[3]
    annot = args[4]

    print()
    build_cell2location_model(sc_input, output_dir, use_gpu, annot)
  
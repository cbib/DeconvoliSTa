# functions.py

import subprocess

def build_cell2location_model(sc_input, params):
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
    cuda_device = params['cuda_device'] if params['gpu'] else "cpu"
    output_dir = params.get('output_dir', '.')
    
    print(f"Building cell2location model with {'GPU' if params['gpu'] else 'CPU'}...")
    
    command = [
        "bash", "-c", f"""
        source activate cell2loc_env &&
        python {params['rootdir']}/scripts/build_model.py \\
            {sc_input} {cuda_device} -a {params['annot']} {sample_id_arg} {epochs} {args} -o {output_dir}
        """
    ]
    print(command)
    subprocess.run(command, check=True)

def fit_cell2location_model(sp_input, sp_input_rds, model, params):
    """
    Fit cell2location model.
    
    Args:
        sp_input (str): Path to spatial input file.
        sp_input_rds (str): Path to RDS spatial input file.
        model (str): Path to the model file.
        params (dict): Parameters for fitting the model.
    """
    output_suffix = sp_input.split('/')[-1]
    output = f"proportions_cell2location_{output_suffix}{params['runID_props']}.preformat"
    epochs = f"-e {params['epoch_fit']}" if params['epoch_fit'] != "default" else ""
    args = params.get('deconv_args', {}).get('cell2location', {}).get('fit', "")
    cuda_device = params['cuda_device'] if params['gpu'] else "cpu"
    output_dir = params.get('output_dir', '.')
    
    print(f"Fitting cell2location model from file {model} with {'GPU' if params['gpu'] else 'CPU'}...")
    print(f"Arguments: {args}")
    
    command = [
        "bash", "-c", f"""
        source activate cell2loc_env &&
        python {params['rootdir']}/scripts/fit_model.py \\
            {sp_input} {model} {cuda_device} {epochs} {args} -o {output_dir} &&
        mv {output_dir}/proportions.tsv {output}
        """
    ]
    print(command)
    subprocess.run(command, check=True)



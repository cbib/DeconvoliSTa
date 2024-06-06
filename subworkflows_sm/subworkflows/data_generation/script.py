
import shutil
synthspot_types_map = {
    'aud': "artificial_uniform_distinct", 
    'add': "artificial_diverse_distinct",
    'auo': "artificial_uniform_overlap", 
    'ado': "artificial_diverse_overlap",
    'adcd': "artificial_dominant_celltype_diverse", 
    'apdcd': "artificial_partially_dominant_celltype_diverse",
    'adrcd': "artificial_dominant_rare_celltype_diverse", 
    'arrcd': "artificial_regional_rare_celltype_diverse",
    'prior': "prior_from_data", 
    'real': "real", 
    'rm': "real_missing_celltypes_visium",
    'arm': "artificial_missing_celltypes_visium", 
    'addm': "artificial_diverse_distinct_missing_celltype_sc",
    'adom': "artificial_diverse_overlap_missing_celltype_sc"
}

synthspot_types_fullnames = list(synthspot_types_map.values())
synthspot_types_flat = synthspot_types_flat = [item for sublist in synthspot_types_map.items() for item in sublist] #All key and values in a list

def create_file_if_not_exists(file_path):
    # Obtenir le répertoire du chemin du fichier
    directory = os.path.dirname(file_path)

    # Créer le répertoire s'il n'existe pas
    if directory and not os.path.exists(directory):
        os.makedirs(directory)

    # Créer le fichier s'il n'existe pas
    if not os.path.exists(file_path):
        with open(file_path, 'w') as f:
            pass  # Ne rien écrire, juste créer le fichier vide

def generate_synthetic_data(sc_input, dataset_type, rep, rootdir, outdir, args=None):
    output_file = f"{os.path.splitext(sc_input)[0]}_{dataset_type}_rep{rep}.rds"
    create_file_if_not_exists(output_file)
    args_str = args if args else ''
    shell_command = (
        f"Rscript {rootdir}/subworkflows/data_generation/generate_synthetic_data.R "
        f"--sc_input {sc_input} --dataset_type {dataset_type} --rep {rep} {args_str}"
    )
    print(shell_command)
    os.system(shell_command)
    # Copy the output file to the output directory
    # output_path = os.path.join(outdir, output_file)
    # shutil.copy(output_file, output_path)

def list_to_dict(flat_list):
    flat_list[0] = flat_list[0].replace("{", "")
    flat_list[-1] = flat_list[-1].replace("}", "")
    # Enlever les virgules à la fin de chaque élément
    flat_list = [item.rstrip(',') for item in flat_list]
    # Créer le dictionnaire à partir des paires clé-valeur
    result_dict = {}
    for i in range(0, len(flat_list), 2):
        key = flat_list[i].strip()[:-1]
        value = flat_list[i + 1].strip()
        result_dict[key] = value

    return result_dict

if __name__ == "__main__":
    import ast
    import os
    import sys  
    # Récupérer tous les arguments de la ligne de commande
    args = sys.argv
    # print(args[1:])
    # config = ast.literal_eval(args[1:])*
    config  = list_to_dict(args[1:])
    # print(config)

    sc_input = config["sc_input"]
    dataset_type = config["dataset_type"]
    reps = config["reps"]
    rootdir = config["rootdir"]


    sc_input_type = 'h5ad' if sc_input.endswith(('h5', 'h5ad')) else 'rds'
    print(f"The synthetic data is of {sc_input_type} format.")

    sc_input_conv = sc_input
    if sc_input_type == 'h5ad':
        sc_input_conv = convert_between_rds_and_h5ad(sc_input_conv)
                
    synthspot_type_input = [synthspot_types_map.get(t, t) for t in config['dataset_type'].split(',') if t in synthspot_types_flat or t in synthspot_types_map]
    synthspot_type_input = list(dict.fromkeys(synthspot_type_input))

    # synthspot_args_input = params.args
    synthspot_args_input =  ' '.join([f"--{k} {v}" for k, v in config.items() if k not in ["dataset_type", "reps", "sc_input"]])

    print(f"Single-cell reference: {sc_input}")
    print(f"Dataset types to be generated: {', '.join(synthspot_type_input)}")
    print(f"Number of replicates per dataset type: {reps}")
    print(f"Arguments: {synthspot_args_input}")

    for dataset_type in synthspot_type_input:
        for rep in range(1, int(reps) + 1):
            output_file = f"{os.path.splitext(os.path.basename(sc_input))[0]}_{dataset_type}_rep{rep}.rds"
            output_path = os.path.join("path/to/synthetic/data", output_file)
            if not os.path.exists(output_path):
                print("this is the fiiiiiile", output_file)
                # Ensure Snakemake knows about the generated files
                shell(f"touch {output_path}")
                generate_synthetic_data(sc_input_conv, dataset_type, rep, rootdir, "path/to/synthetic/data", synthspot_args_input)

import os
import sys
import yaml

# Lire le fichier de configuration YAML
with open("config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

# Fonction pour obtenir le nom de base du fichier sans extension
def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

# Helper function to check if a config variable is present
def get_config_var(config, var_name, default=None):
    if var_name not in config:
        if default is None:
            raise ValueError(f"Error: '{var_name}' is missing in the configuration file.")
        else:
            print(f"Warning: '{var_name}' is missing in the configuration file. Using default: {default}")
            return default
    return config[var_name]

# Charger les paramètres de configuration nécessaires

mode = get_config_var(config, "mode")
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

skip_metrics = get_config_var(config, "skip_metrics", "false")

if mode == "run_dataset":
    sc_input = get_config_var(config, "sc_input")
    sp_input = get_config_var(config, "sp_input")
    # methods = get_config_var(config, "methods").split(',')
    supported_methods = {"cell2location", "rctd", "spatialdwls", "nnls", "ddls", "dirichlet", "seurat"}
    methods = get_config_var(config, "methods").split(',')
    methods_set = set(methods)
    # Find unsupported methods by subtracting the supported methods from the methods set
    unsupported_methods = methods_set - supported_methods
    # If there are unsupported methods, raise an error with the appropriate message
    if unsupported_methods:
        raise ValueError(f"Error: The following methods are not supported: {', '.join(unsupported_methods)}")

    output_suffix = get_basename(sp_input)
    runID_props = get_config_var(params, "runID_props")
    include: "subworkflows/deconvolution/run_methods.smk"
    output_dir = get_config_var(config, "output")

    output_files = [f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv" for method in methods]
    if skip_metrics == "false":
        include: "subworkflows/evaluation/evaluate_methods.smk"
        metrics_files = [f"{output_dir}/metrics/metrics_{method}_{output_suffix}{runID_props}.tsv" for method in methods]
        rule main:
            input:
                output_files,
                metrics_files
    else:
        rule main:
            input:
                output_files
elif mode == "generate_data":
    sc_input = get_config_var(config, "sc_input")
    dataset_types = [synthspot_types_map[t] for t in get_config_var(config, "dataset_type").split(',')]
    reps = get_config_var(config, "reps")
    rootdir = get_config_var(config, "rootdir")

    generated_files = [f"synthetic_data_sm/{os.path.basename(sc_input).split('.')[0]}_{dataset_type}_rep{rep}.rds" for dataset_type in dataset_types for rep in range(1, int(reps) + 1)]
    include: "subworkflows/data_generation/generate_data.smk"
    rule gen_files:
        input:
            generated_files
elif mode == "generate_vis":
    sp_input = get_config_var(config, "sp_input")
    output_dir = get_config_var(config, "output", ".")
    generated_file = f"{output_dir}/{os.path.basename(sp_input).split('.')[0]}.html"
    
    norm_weights_filepaths = get_config_var(config, "norm_weights_filepaths").split(",")
    raw_norm_weights_filepaths_without_split = get_config_var(config, "norm_weights_filepaths")
    st_coords_filepath = get_config_var(config, "st_coords_filepath")
    data_clustered = get_config_var(config, "data_clustered")
    image_path = get_config_var(config, "image_path")
    n_largest_cell_types = get_config_var(config, "n_largest_cell_types", "5")
    scale_factor = get_config_var(config, "scale_factor")
    deconv_methods = get_config_var(config, "deconv_methods")

    print(raw_norm_weights_filepaths_without_split)
    rule gen_html:
        input:
            norm_weights_filepaths = norm_weights_filepaths,
            st_coords_filepath = st_coords_filepath,
            data_clustered = data_clustered,
            image_path = image_path
        output:
            generated_file = generated_file
        shell:
            """
            start_time=$(date +%s)
            python3 subworkflows/visualization/sp_visualizer.py {sp_input} {raw_norm_weights_filepaths_without_split} {st_coords_filepath} {data_clustered} {image_path} {n_largest_cell_types} {scale_factor} {generated_file} {deconv_methods}
            end_time=$(date +%s)
            elapsed_time=$((end_time - start_time))
            echo "html_generation took $elapsed_time seconds"
            """
else:
    raise ValueError("Error: Enter a valid execution mode --mode")

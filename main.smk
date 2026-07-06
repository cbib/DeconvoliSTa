import os
import sys
import yaml

# Read the YAML configuration file
with open("config.yaml", "r") as config_file:
    params = yaml.safe_load(config_file)

# Function to get the file basename (without extension)
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

# Load the required configuration parameters

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
    methods = [m.strip() for m in get_config_var(config, "methods").split(',') if m.strip()]
    # A method is supported iff its subworkflow exists (subworkflows/deconvolution/<m>/run_<m>.smk).
    # This auto-discovers methods: adding one needs no edit here (matches ./deconvolista discovery).
    unsupported_methods = [m for m in methods
                           if not os.path.exists(f"subworkflows/deconvolution/{m}/run_{m}.smk")]
    if unsupported_methods:
        raise ValueError(
            "Error: unknown method(s) — no subworkflows/deconvolution/<m>/run_<m>.smk for: "
            f"{', '.join(unsupported_methods)}")

    output_suffix = get_basename(sp_input)
    runID_props = get_config_var(params, "runID_props")
    include: "subworkflows/deconvolution/run_methods.smk"
    output_dir = get_config_var(config, "output")
    do_visu = get_config_var(config, "do_visu", "false")

    output_files = [f"{output_dir}/proportions_{method}_{output_suffix}{runID_props}.tsv" for method in methods]
    main_inputs = list(output_files)

    if skip_metrics == "false":
        include: "subworkflows/evaluation/evaluate_methods.smk"
        metrics_files = [f"{output_dir}/metrics/metrics_{method}_{output_suffix}{runID_props}.tsv" for method in methods]
        main_inputs += metrics_files

    # do_visu=true: from the spatial Seurat object, chain BayesSpace clustering + tissue
    # image/coordinate extraction + the interactive HTML (no manually-prepared inputs).
    if do_visu == "true":
        include: "subworkflows/visualization/auto_visualize.smk"
        main_inputs += [_html]

    rule main:
        input:
            main_inputs
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

import os

def get_config_var(config, var_name, default=None):
    if var_name not in config:
        if default is None:
            raise ValueError(f"Missing '{var_name}'")
        print(f"Using default for {var_name}: {default}")
        return default
    return config[var_name]

sp_input = get_config_var(config, "sp_input")
output_dir = get_config_var(config, "output", ".")
generated_file = f"{output_dir}/{os.path.basename(sp_input).split('.')[0]}.html"

norm_weights = get_config_var(config, "norm_weights_filepaths").split(",")
raw_norm_weights = get_config_var(config, "norm_weights_filepaths")
st_coords = get_config_var(config, "st_coords_filepath")
data_clustered = get_config_var(config, "data_clustered")
image_path = get_config_var(config, "image_path")
n_largest_cell_types = get_config_var(config, "n_largest_cell_types", "5")
scale_factor = get_config_var(config, "scale_factor")
deconv_methods = get_config_var(config, "deconv_methods")

rule all:
    input:
        generated_file

rule generate_html:
    input:
        norm_weights_filepaths=norm_weights,
        st_coords_filepath=st_coords,
        data_clustered=data_clustered,
        image_path=image_path
    output:
        html=generated_file
    params:
        sp_input=sp_input,
        raw_norm_weights=raw_norm_weights,
        st_coords=st_coords,
        data_clustered=data_clustered,
        image_path=image_path,
        n_largest_cell_types=n_largest_cell_types,
        scale_factor=scale_factor,
        deconv_methods=deconv_methods
    shell:
        """
        start=$(date +%s)
        python3 subworkflows/visualization/sp_visualizer.py \
            "{params.sp_input}" "{params.raw_norm_weights}" "{params.st_coords}" \
            "{params.data_clustered}" "{params.image_path}" "{params.n_largest_cell_types}" \
            "{params.scale_factor}" "{output.html}" "{params.deconv_methods}"
        end=$(date +%s)
        echo "html_generation took $((end-start)) seconds"
        """

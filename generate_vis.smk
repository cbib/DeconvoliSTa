import os, sys, yaml

def get_config_var(config, var_name, default=None):
    if var_name not in config:
        if default is None:
            raise ValueError(f"Missing '{var_name}'")
        print(f"Using default for {var_name}: {default}")
        return default
    return config[var_name]

with open("config.yaml", "r") as f:
    cfg = yaml.safe_load(f)

sp_input          = get_config_var(cfg, "sp_input")
output_dir        = get_config_var(cfg, "output", ".")
generated_file    = f"{output_dir}/{os.path.basename(sp_input).split('.')[0]}.html"

norm_weights      = get_config_var(cfg, "norm_weights_filepaths").split(",")
raw_norm_weights  = get_config_var(cfg, "norm_weights_filepaths")   # keep original string
st_coords         = get_config_var(cfg, "st_coords_filepath")
data_clustered    = get_config_var(cfg, "data_clustered")
image_path        = get_config_var(cfg, "image_path")
n_largest_cell_types = get_config_var(cfg, "n_largest_cell_types", "5")
scale_factor      = get_config_var(cfg, "scale_factor")
deconv_methods    = get_config_var(cfg, "deconv_methods")

rule all:
    input: generated_file

rule generate_html:
    input:
        norm_weights_filepaths=norm_weights,
        st_coords_filepath=st_coords,
        data_clustered=data_clustered,
        image_path=image_path
    output:
        html=generated_file
    shell:
        """
        start=$(date +%s)
        python3 subworkflows/visualization/sp_visualizer.py \
            {sp_input} {raw_norm_weights} {st_coords_filepath} \
            {data_clustered} {image_path} {n_largest_cell_types} \
            {scale_factor} {generated_file} {deconv_methods}
        end=$(date +%s)
        echo "html_generation took $((end-start)) seconds"
        """
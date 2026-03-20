# subworkflows/visualization/sp_visualizer.smk
# This Snakefile launches the Python visualization script `sp_visualizer.py`.
# It is designed to be included from a higher‑level workflow (e.g. generate_vis.smk).
# The rule simply passes through all required inputs and executes the python command.

import os

rule sp_visualize:
    input:
        # Path to the spatial data file
        sp_input=os.path.expanduser("{sp_input}"),
        # Comma separated list of norm weight files
        norm_weights_filepaths=os.path.expanduser("{norm_weights_filepaths}"),
        st_coords_filepath=os.path.expanduser("{st_coords_filepath}"),
        data_clustered=os.path.expanduser("{data_clustered}"),
        image_path=os.path.expanduser("{image_path}")
    output:
        html=os.path.expanduser("{generated_file}")
    params:
        n_largest_cell_types="{n_largest_cell_types}",
        scale_factor="{scale_factor}",
        deconv_methods="{deconv_methods}"
    shell:
        """
        start=$(date +%s)
        python3 subworkflows/visualization/sp_visualizer.py \
            {sp_input} {norm_weights_filepaths} {st_coords_filepath} \
            {data_clustered} {image_path} {params.n_largest_cell_types} \
            {params.scale_factor} {output.html} {params.deconv_methods}
        end=$(date +%s)
        echo "html_generation took $((end-start)) seconds"
        """
# subworkflows/visualization/sp_visualizer.smk
# This Snakefile launches the Python visualization script `sp_visualizer.py

# Rule: generate_vis

## Purpose  
Builds an interactive HTML report overlaying deconvolution proportions on a Visium image and displaying spot clustering.

### Required configuration
| Key                        | Type    | Example                                  |
|----------------------------|---------|------------------------------------------|
| mode                       | string  | `"generate_vis"`                         |
| sp_input                   | file    | `"sample_sp.rds"`                        |
| norm_weights_filepaths     | csv list| `"props_rctd.tsv,props_cell2location.tsv"` |
| st_coords_filepath         | CSV file| `"tissue_positions_list.csv"`            |
| data_clustered             | CSV file| `"seurat_metadata.csv"`                 |
| image_path                 | image   | `"tissue_hires.png"`                    |
| deconv_methods             | csv list| `"rctd,cell2location"`                  |

### Optional configuration
| Key                     | Default | Notes                                            |
|-------------------------|---------|--------------------------------------------------|
| output                  | `.`     | Directory where the HTML file will be written.   |
| n_largest_cell_types    | `"5"`   | How many top cell types to show per spot.       |
| scale_factor            | `1.0`   | Scaling factor for pixel coordinates.           |

### Example
```bash
snakemake -s generate_vis.smk \
  --config mode="generate_vis" sp_input="UKF243_T_ST_1_raw.rds" output="vis_output" \
  norm_weights_filepaths="props_rctd.tsv,props_cell2location.tsv" \
  st_coords_filepath="tissue_positions_list_243.csv" data_clustered="seurat_metadata.csv" \
  image_path="tissue_hires.png" scale_factor=0.24414062 deconv_methods=rctd,cell2location
```

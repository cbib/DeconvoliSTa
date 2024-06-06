# Snakefile

import os
import sys
sys.path.append('cell2location/')
from functions import build_cell2location_model, fit_cell2location_model


# Préparez les chemins d'entrée/sortie
sc_input = config["sc_input"]
sp_input = config["sp_input"]


# rule all:
#     input:
#         "sc.h5ad",
#         # "proportions_cell2location_{output_suffix}{runID_props}.preformat"

rule build_cell2location:
    input:
        sc_input
    output:
        "sc.h5ad"
    run:
        build_cell2location_model(input.sc_input, config)

rule fit_cell2location:
    input:
        sp_input,
        model="sc.h5ad"
    output:
        "proportions_cell2location_{output_suffix}{runID_props}.preformat"
    run:
        fit_cell2location_model(input.sp_input, input.sp_input, input.model, config)

rule format_c2l:
    input:
        "proportions_cell2location_{output_suffix}{runID_props}.preformat"
    output:
        "formatted_proportions_cell2location_{output_suffix}{runID_props}.tsv"
    run:
        format_tsv(input.input, output.output, config)




# Snakefile

import os
import sys
from functions import build_cell2location_model, fit_cell2location_model
# sys.path.append('/subworkflows_sm/subworkflows/data_generation/')

# configfile: "config.yaml"

# Préparez les chemins d'entrée/sortie
sc_input = config["sc_input"]
sp_input = config["sp_input"]
sp_input_rds = config["sp_input_rds"]

rule all:
    input:
        "sc.h5ad",
        "proportions_cell2location_{output_suffix}{runID_props}.preformat"

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
        sp_input_rds,
        model="sc.h5ad"
    output:
        "proportions_cell2location_{output_suffix}{runID_props}.preformat"
    run:
        fit_cell2location_model(input.sp_input, input.sp_input_rds, input.model, config)

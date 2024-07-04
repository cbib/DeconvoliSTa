# Snakefile
import os
import yaml


import sys
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

# Print the full Python version
print(f"Full Python version: {sys.version}")
rule generateSyntheticData:
    input:
        sc_input=config['sc_input']
    output:
        expand("synthetic_data_sm/{basename}_{dataset_type}_rep{rep}.rds", 
               basename= os.path.splitext(os.path.basename(config['sc_input']))[0],
               dataset_type = [synthspot_types_map[t] for t in config['dataset_type'].split(',')],
               rep=range(1, int(config['reps']) + 1))
    singularity: 
        "docker://csangara/synthspot:latest"
    params:
        rootdir=config['rootdir'],
        reps=config['reps'],
        args= ' '.join([f"--{k} {v}" for k, v in config.items() if k not in ["dataset_type", "reps", "sc_input"]])
    shell:
        "python3  ./subworkflows_sm/data_generation/script.py {config}"
            

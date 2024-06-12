# Snakefile
import os

import sys

# Print the full Python version
print(f"Full Python version: {sys.version}")

rule generateSyntheticData:
    input:
        sc_input=config['sc_input']
    output:
        expand("synthetic_data_sm/{basename}_{dataset_type}_rep{rep}.rds", 
               basename= os.path.splitext(os.path.basename(config['sc_input']))[0],
               dataset_type = config['dataset_type'].split(','),
               rep=range(1, int(config['reps']) + 1))
    singularity: 
        "docker://csangara/synthspot:latest"
    params:
        rootdir=config['rootdir'],
        reps=config['reps'],
        args= ' '.join([f"--{k} {v}" for k, v in config.items() if k not in ["dataset_type", "reps", "sc_input"]])
    shell:
        "python3  ./subworkflows_sm/data_generation/script.py {config}"
            

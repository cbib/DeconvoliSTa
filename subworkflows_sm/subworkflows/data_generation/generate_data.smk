# include: "../helper_processes.smk"
# rule convert_sc = rules.convert_between_rds_and_h5ad

# Initialization: possible dataset types
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

synthspot_types_fullnames = list(synthspot_types_map.values())
synthspot_types_flat = synthspot_types_flat = [item for sublist in synthspot_types_map.items() for item in sublist] #All key and values in a list

 
import os
import shutil

def generate_synthetic_data(sc_input, dataset_type, rep, rootdir, outdir, args=None):
    output_file = f"{os.path.splitext(sc_input)[0]}_{dataset_type}_rep{rep}.rds"
    args_str = args if args else ''
    shell_command = (
        f"Rscript {rootdir}/subworkflows/data_generation/generate_synthetic_data.R "
        f"--sc_input {sc_input} --dataset_type {dataset_type} --rep {rep} {args_str}"
    )
    print(shell_command)
    os.system(shell_command)
    # Copy the output file to the output directory
    output_path = os.path.join(outdir, output_file)
    shutil.copy(output_file, output_path)

rule generateSyntheticData:
    input:
        sc_input=config['sc_input']
    output:
        expand("synthetic/data/{basename}_{dataset_type}_rep{rep}.rds", 
               basename= os.path.splitext(os.path.basename(config['sc_input']))[0],
               dataset_type=config['dataset_type'].split(','), 
               rep=range(1, int(config['reps']) + 1))
    params:
        rootdir=config['rootdir'],
        reps=config['reps'],
        args= ' '.join([f"--{k} {v}" for k, v in config.items() if k not in ["dataset_type", "reps", "sc_input"]])
    container: 
        "csangara/synthspot:latest"
    shell:
        "python3  ./subworkflows_sm/subworkflows/data_generation/script.py {config}"
            

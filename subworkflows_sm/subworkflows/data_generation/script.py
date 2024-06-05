# def main():
#     return 
if __name__ == "__main__":
    import ast

    import sys

    # Récupérer tous les arguments de la ligne de commande
    args = sys.argv
    print(args[1:])
    config = ast.literal_eval(args[1:])

    sc_input = config["sc_input"]
    dataset_type = config["dataset_type"]
    rep = config["reps"]
    rootdir = config["rootdir"]


    sc_input_type = 'h5ad' if sc_input.endswith(('h5', 'h5ad')) else 'rds'
    print(f"The synthetic data is of {sc_input_type} format.")

    sc_input_conv = sc_input
    if sc_input_type == 'h5ad':
        sc_input_conv = convert_between_rds_and_h5ad(sc_input_conv)
                
    synthspot_type_input = [synthspot_types_map.get(t, t) for t in config['dataset_type'].split(',') if t in synthspot_types_flat or t in synthspot_types_map]
    synthspot_type_input = list(dict.fromkeys(synthspot_type_input))

    # synthspot_args_input = params.args
    synthspot_args_input =  params.args # a changer

    print(f"Single-cell reference: {sc_inputt}")
    print(f"Dataset types to be generated: {', '.join(synthspot_type_input)}")
    print(f"Number of replicates per dataset type: {params.reps}")
    print(f"Arguments: {synthspot_args_input}")

    for dataset_type in synthspot_type_input:
        for rep in range(1, int(params.reps) + 1):
            output_file = f"{os.path.splitext(os.path.basename(sc_input))[0]}_{dataset_type}_rep{rep}.rds"
            output_path = os.path.join("path/to/synthetic/data", output_file)
            if not os.path.exists(output_path):
                generate_synthetic_data(sc_input_conv, dataset_type, rep, rootdir, "path/to/synthetic/data", synthspot_args_input)
                # Ensure Snakemake knows about the generated files
                print(output_file)
                shell(f"touch {output_path}")

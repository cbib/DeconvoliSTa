
include: "../helper_processes.smk"
rule convert_sc = rules.convert_between_rds_and_h5ad

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

rule test: 
    print(synthspot_types_fullnames)
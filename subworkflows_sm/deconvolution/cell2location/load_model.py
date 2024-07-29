#!/usr/bin/env python3
import argparse as arp
import os

def main():
    ##### PARSING COMMAND LINE ARGUMENTS #####
    prs = arp.ArgumentParser()
    
    prs.add_argument('sp_data_path', type = str, help = 'path of spatial data')
    prs.add_argument('model_path', type = str, help = 'path to trained cell2location model')
    prs.add_argument('cuda_device', type = str, help = "index of cuda device ID or cpu")
    prs.add_argument('-o','--out_dir', default = os.getcwd(), type = str, help = 'proportions output directory')
    prs.add_argument('-p', '--posterior_sampling', default=1000, type = int, help = "number of samples to take from the posterior distribution")
    
    args = prs.parse_args()
    
    cuda_device = args.cuda_device
    print("cuda device = ", cuda_device)

    sp_data_path = args.sp_data_path
    output_folder = args.out_dir
        
    assert (cuda_device.isdigit() or cuda_device == "cpu"), "invalid device input"
    
    ##### MAIN PART #####
    if cuda_device.isdigit():
        os.environ["CUDA_VISIBLE_DEVICES"]=cuda_device

    import scanpy as sc
    import cell2location
    import scvi
    import anndata

    # silence scanpy that prints a lot of warnings
    import warnings
    warnings.filterwarnings('ignore')

    print("Reading in spatial data from " + sp_data_path + "...")
    adata = sc.read_h5ad(sp_data_path)
    adata.var['SYMBOL'] = adata.var_names

    # mitochondria-encoded (MT) genes should be removed for spatial mapping
    adata.var['mt'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]
    adata = adata[:, ~adata.var['mt'].values].copy()

    # # prepare anndata for cell2location model
    # scvi.data.setup_anndata(adata=adata)
    print("printing model")
    print(adata.obs)
    print("printing model")
    sc_data = anndata.read_h5ad('/home/abderahim/spotless-benchmark/test_sc_data.h5ad')
    print(sc_data.obs["subclass"])
    adata.obs["subclass"] =  ['Astro', 'CR', 'Endo', 'L2/3 IT', 'L4', 'L5 IT', 'L5 PT', 'L6 CT', 'L6 IT', 'L6b' ,\
        'Lamp5', 'Macrophage', 'Meis2', 'NP' ,'Oligo', 'Peri' ]
    scvi.data.setup_anndata(
        adata=adata,
        labels_key="subclass", 
        batch_key=None,  
        categorical_covariate_keys=None, 
    )

    # Load the trained model
    mod = cell2location.models.Cell2location.load("/home/abderahim/spotless-benchmark/res", adata=adata)

    # Export the estimated cell abundance (summary of the posterior distribution)
    adata = mod.export_posterior(
        adata, sample_kwargs={'num_samples': args.posterior_sampling,
        'batch_size': mod.adata.n_obs, 'use_gpu': cuda_device.isdigit()}
    )

    # Export proportion file, but first rename columns and divide by rowSums
    props = adata.obsm['q05_cell_abundance_w_sf']
    props = props.rename(columns={x:x.replace("q05cell_abundance_w_sf_", "") for x in props.columns})
    props = props.div(props.sum(axis=1), axis='index')
    props.to_csv(os.path.join(output_folder, 'proportions.tsv'), sep="\t")

if __name__ == '__main__':
    main()

# #!/usr/bin/env python3
# import argparse as arp
# import os

# def main():
#     ##### PARSING COMMAND LINE ARGUMENTS #####
#     prs = arp.ArgumentParser()
    
#     prs.add_argument('sp_data_path', type=str, help='path of spatial data')
#     prs.add_argument('model_path', type=str, help='path to trained model (.h5ad)')
#     prs.add_argument('cuda_device', type=str, help='index of cuda device ID or cpu')
#     prs.add_argument('-o', '--out_dir', default=os.getcwd(), type=str, help='output directory for results')

#     args = prs.parse_args()
    
#     cuda_device = args.cuda_device
#     print("cuda device =", cuda_device)

#     sp_data_path = args.sp_data_path
#     output_folder = args.out_dir
        
#     assert (cuda_device.isdigit() or cuda_device == "cpu"), "invalid device input"
    
#     if cuda_device.isdigit():
#         os.environ["CUDA_VISIBLE_DEVICES"] = cuda_device

#     import scanpy as sc
#     import anndata
#     import pandas as pd
#     import numpy as np

#     import cell2location
#     import scvi

#     import matplotlib as mpl
#     from matplotlib import rcParams
#     rcParams['pdf.fonttype'] = 42
#     import matplotlib.pyplot as plt
#     import seaborn as sns

#     # silence scanpy that prints a lot of warnings
#     import warnings
#     warnings.filterwarnings('ignore')

#     print("Reading in spatial data from " + sp_data_path + "...")
#     adata = sc.read_h5ad(sp_data_path)
#     adata.var['SYMBOL'] = adata.var_names

#     # mitochondria-encoded (MT) genes should be removed for spatial mapping
#     adata.var['mt'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]
#     adata = adata[:, ~adata.var['mt'].values]

#     adata_vis = adata.copy()
#     adata_vis.raw = adata_vis

#     print("Reading in the model...")
#     adata_scrna_raw = sc.read("/home/abderahim/spotless-benchmark/res/sc.h5ad")
    
#     # Export estimated expression in each cluster
#     if 'means_per_cluster_mu_fg' in adata_scrna_raw.varm.keys():
#         inf_aver = adata_scrna_raw.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
#                                         for i in adata_scrna_raw.uns['mod']['factor_names']]].copy()
#     else:
#         inf_aver = adata_scrna_raw.var[[f'means_per_cluster_mu_fg_{i}' 
#                                         for i in adata_scrna_raw.uns['mod']['factor_names']]].copy()
#     inf_aver.columns = adata_scrna_raw.uns['mod']['factor_names']

#     # find shared genes and subset both anndata and reference signatures
#     intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
#     adata_vis = adata_vis[:, intersect].copy()
#     inf_aver = inf_aver.loc[intersect, :].copy()

#     print("Loading the trained model...")
#     trained_model = cell2location.models.Cell2location.load("/home/abderahim/spotless-benchmark/res", adata_vis)
#     # prepare anndata for cell2location model
#     scvi.data.setup_anndata(adata=adata_vis)

#     print("Predicting cell abundance...")
#     adata_vis = trained_model.export_posterior(
#         adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': adata_vis.n_obs, 'use_gpu': cuda_device.isdigit()}
#     )

#     # Save anndata object with results
#     adata_vis.write(os.path.join(output_folder, 'sp.h5ad'))

#     # Export proportion file, but first rename columns and divide by rowSums
#     props = adata_vis.obsm['q05_cell_abundance_w_sf']
#     props = props.rename(columns={x: x.replace("q05cell_abundance_w_sf_", "") for x in props.columns})
#     props = props.div(props.sum(axis=1), axis='index')
#     props.to_csv(os.path.join(output_folder, 'proportions.tsv'), sep="\t")

# if __name__ == '__main__':
#     main()

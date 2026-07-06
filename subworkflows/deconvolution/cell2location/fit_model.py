#!/usr/bin/env python3
import argparse as arp
import os


try:
    from pybiomart import Server
except ImportError:
    Server = None


def convert_query_geneSymbol_to_ensemblID(adata):
    import pandas as pd

    if Server is None:
        raise ImportError("pybiomart is required only when map_genes=true")

    server = Server(host='http://www.ensembl.org')
    dataset = (
        server.marts['ENSEMBL_MART_ENSEMBL']
        .datasets['hsapiens_gene_ensembl']
    )

    gene_map = dataset.query(
        attributes=['ensembl_gene_id', 'external_gene_name']
    )

    adata.var['SYMBOL'] = adata.var_names
    adata = adata[:, adata.var['SYMBOL'].isin(gene_map['Gene name'].values)].copy()
    adata.var['ENSEMBL'] = adata.var['SYMBOL'].map(
        dict(zip(gene_map['Gene name'], gene_map['Gene stable ID']))
    )
    adata.var_names = adata.var['ENSEMBL']
    adata.var['features'] = adata.var['ENSEMBL']
    return adata


def main():
    prs = arp.ArgumentParser()

    prs.add_argument('sp_data_path', type=str, help='path of spatial data')
    prs.add_argument('model_path', type=str, help='path to regression model')
    prs.add_argument('cuda_device', type=str, help='index of cuda device ID or cpu')

    prs.add_argument(
        '-o', '--out_dir',
        default=os.getcwd(),
        type=str,
        help='model and proportions output directory'
    )
    prs.add_argument(
        '-e', '--epochs',
        default=30000,
        type=int,
        help='number of epochs to fit the model'
    )
    prs.add_argument(
        '-p', '--posterior_sampling',
        default=1000,
        type=int,
        help='number of samples to take from the posterior distribution'
    )
    prs.add_argument(
        '-n', '--n_cells_per_location',
        default=8,
        type=int,
        help='estimated number of cells per spot'
    )
    prs.add_argument(
        '-d', '--detection_alpha',
        default=200,
        type=int,
        help='within-experiment variation in RNA detection sensitivity'
    )
    prs.add_argument(
        '-m', '--map_genes',
        default="false",
        type=str,
        help='map genes between single cell and spatial'
    )

    args = prs.parse_args()

    cuda_device = args.cuda_device
    print("cuda device requested =", cuda_device)

    sp_data_path = args.sp_data_path
    output_folder = args.out_dir

    assert (cuda_device.isdigit() or cuda_device == "cpu"), "invalid device input"
    accelerator = "cpu" if cuda_device == "cpu" else "gpu"
    print("Using accelerator:", accelerator)

    print("Parameters\n==========")
    print("Detection alpha: {}\nCells per location: {}".format(
        args.detection_alpha, args.n_cells_per_location
    ))
    print("Training epochs: {}\nPosterior sampling: {}".format(
        args.epochs, args.posterior_sampling
    ))
    print("Map genes: {}".format(args.map_genes))
    print("==========")

    import scanpy as sc
    import numpy as np
    import cell2location
    from cell2location.models import Cell2location

    import warnings
    warnings.filterwarnings('ignore')

    print("Reading in spatial data from " + sp_data_path + "...")
    adata = sc.read_h5ad(sp_data_path)
    adata.var['SYMBOL'] = adata.var_names

    adata.var['mt'] = [gene.lower().startswith('mt-') for gene in adata.var['SYMBOL']]
    adata = adata[:, ~adata.var['mt'].values].copy()

    adata_vis = adata.copy()
    adata_vis.raw = adata_vis

    if str(args.map_genes).lower() == "true":
        adata_vis = convert_query_geneSymbol_to_ensemblID(adata_vis)

    print("Reading in the model...")
    adata_scrna_raw = sc.read(args.model_path)

    if 'means_per_cluster_mu_fg' in adata_scrna_raw.varm.keys():
        inf_aver = adata_scrna_raw.varm['means_per_cluster_mu_fg'][
            [f'means_per_cluster_mu_fg_{i}' for i in adata_scrna_raw.uns['mod']['factor_names']]
        ].copy()
    else:
        inf_aver = adata_scrna_raw.var[
            [f'means_per_cluster_mu_fg_{i}' for i in adata_scrna_raw.uns['mod']['factor_names']]
        ].copy()

    inf_aver.columns = adata_scrna_raw.uns['mod']['factor_names']

    print("Finding shared genes...")
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    print("Preparing anndata for cell2location...")
    Cell2location.setup_anndata(adata=adata_vis)

    print("Building cell2location model...")
    mod = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=args.n_cells_per_location,
        detection_alpha=args.detection_alpha
    )

    print("Training cell2location model on", accelerator, "...")
    mod.train(
        max_epochs=args.epochs,
        batch_size=None,
        train_size=1,
        accelerator=accelerator
    )

    print("Exporting posterior on CPU...")
    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={
            'num_samples': args.posterior_sampling,
            'batch_size': mod.adata.n_obs
        }
    )

    print("Writing proportions.tsv ...")
    props = adata_vis.obsm['q05_cell_abundance_w_sf'].copy()
    props = props.rename(
        columns={x: x.replace("q05cell_abundance_w_sf_", "") for x in props.columns}
    )
    props = props.div(props.sum(axis=1), axis='index')

    props.to_csv(
        os.path.join(output_folder, 'proportions.tsv'),
        sep="\t",
        index=True
    )

    print("Done.")


if __name__ == '__main__':
    main()
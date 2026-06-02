#!/usr/bin/env python3

import argparse as arp
import os
import matplotlib as mpl

mpl.use('Agg')


def get_basename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]


def main():
    prs = arp.ArgumentParser()

    prs.add_argument(
        'sc_data_path',
        type=str,
        help='path to single cell h5ad count data'
    )
    prs.add_argument(
        'sp_data_path',
        type=str,
        help='path to spatial data'
    )
    prs.add_argument(
        'cuda_device',
        type=str,
        help='index of cuda device ID or cpu'
    )
    prs.add_argument(
        '-o', '--out_dir',
        default=None,
        type=str,
        help='directory for regression model'
    )
    prs.add_argument(
        '-a', '--annotation_column',
        default='celltype',
        type=str,
        help='column name for covariate'
    )
    prs.add_argument(
        '-s', '--sample_column',
        default=None,
        type=str,
        help='column containing sample id or batch information'
    )
    prs.add_argument(
        '-t', '--tech_column',
        default=None,
        nargs='+',
        type=str,
        help='multiplicative technical effects, such as platform effects'
    )
    prs.add_argument(
        '-e', '--epochs',
        default=250,
        type=int,
        help='number of epochs to train the model'
    )
    prs.add_argument(
        '-p', '--posterior_sampling',
        default=1000,
        type=int,
        help='number of samples to take from the posterior distribution'
    )

    args = prs.parse_args()
    cuda_device = args.cuda_device
    print("cuda device requested = ", cuda_device)
    print("Forcing CPU execution for compatibility with cluster driver.")

    assert (cuda_device.isdigit() or cuda_device == "cpu"), "invalid device id"
    assert os.path.exists(args.sc_data_path), f"{args.sc_data_path} sc file not found"

    if args.out_dir is None:
        output_folder = os.path.dirname(args.sc_data_path) + '/cell2location_results/'
    else:
        output_folder = args.out_dir

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    print("Parameters\n==========")
    print("Training epochs: {}\nPosterior sampling: {}".format(args.epochs, args.posterior_sampling))
    print("==========")

    import anndata
    import scanpy as sc  # noqa: F401
    import pandas as pd  # noqa: F401
    import numpy as np  # noqa: F401

    import cell2location  # noqa: F401
    import scvi  # noqa: F401

    import matplotlib.pyplot as plt  # noqa: F401
    from matplotlib import rcParams
    rcParams['pdf.fonttype'] = 42
    import seaborn as sns  # noqa: F401

    from cell2location.utils.filtering import filter_genes
    from cell2location.models import RegressionModel

    import warnings
    warnings.filterwarnings('ignore')

    print("Reading scRNA-seq data from " + args.sc_data_path + "...")
    adata_scrna_raw = anndata.read_h5ad(args.sc_data_path)
    adata_scrna_raw.var['SYMBOL'] = adata_scrna_raw.var_names

    print("Before filtering: {} cells and {} genes.".format(*adata_scrna_raw.shape))
    selected = filter_genes(
        adata_scrna_raw,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12
    )
    print("After selecting genes.")
    adata_scrna_raw = adata_scrna_raw[:, selected].copy()

    print("After filtering: {} cells and {} genes.".format(*adata_scrna_raw.shape))
    print("Preparing anndata for the regression model...")

    RegressionModel.setup_anndata(
        adata=adata_scrna_raw,
        batch_key=args.sample_column,
        labels_key=args.annotation_column
    )

    print("Building regression model...")
    mod = RegressionModel(adata_scrna_raw)

    print("Training regression model on CPU...")
    mod.train(
        max_epochs=args.epochs,
        batch_size=2500,
        train_size=1,
        lr=0.002,
        accelerator="cpu"
    )

    print("Exporting posterior on CPU...")
    adata_scrna_raw = mod.export_posterior(
        adata_scrna_raw,
        sample_kwargs={
            'num_samples': args.posterior_sampling,
            'batch_size': 2500
        }
    )

    print("Saving model...")
    mod.save(output_folder, overwrite=True)

    output_file = output_folder + f"/sc_{get_basename(args.sc_data_path)}_{get_basename(args.sp_data_path)}.h5ad"

    try:
        adata_scrna_raw.write(output_file)
    except ValueError:
        print("There seems to be an issue with the conversion. Renaming columns...")
        if os.path.exists(output_file):
            os.remove(output_file)
        adata_scrna_raw.__dict__['_raw'].__dict__['_var'] = (
            adata_scrna_raw.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
        )
        adata_scrna_raw.write(output_file)

    print("Done.")


if __name__ == '__main__':
    main()
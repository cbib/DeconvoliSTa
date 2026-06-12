# Automatic visualization chain (enabled with do_visu=true in run_dataset mode).
#
# From the deconvolution proportions + the spatial Seurat object, build the interactive HTML
# with NO manually-prepared inputs:
#   1. extract the tissue image, spot coordinates and scale factor from the Seurat object,
#   2. compute a BayesSpace spatial clustering,
#   3. generate the interactive visualization.
#
# Included from main.smk, which provides: sp_input, output_dir, output_suffix, methods,
# runID_props. The clustering + visualization images are custom (not on a public registry);
# build them once via subworkflows/clustering/bayesspace.def and subworkflows/visualization/visu.def
# (see SETUP.md). do_visu is therefore opt-in; the deconvolution itself needs no custom image.

_vis_dir    = f"{output_dir}/visualization"
_assets_dir = f"{_vis_dir}/assets"
_clustering = f"{_vis_dir}/clustering_bayesspace.csv"   # spatial domains (BayesSpace)
_seurat_cl  = f"{_vis_dir}/clustering_seurat.csv"       # transcriptomic clusters (Seurat)
_html       = f"{_vis_dir}/{output_suffix}_visualization.html"

_bayes_q   = get_config_var(config, "bayes_q", "auto")
_seurat_res = get_config_var(config, "seurat_res", "0.5")   # lower -> fewer, broader clusters
_n_largest = get_config_var(config, "n_largest_cell_types", "5")

_bayes_sif = f"{config.get('sif_dir', 'sif')}/sp_bayesspace.sif"
_visu_sif  = f"{config.get('sif_dir', 'sif')}/visu.sif"

# Proportions in the SAME order as `methods`, so the visualization columns line up with the methods.
_props = [f"{output_dir}/proportions_{m}_{output_suffix}{runID_props}.tsv" for m in methods]

rule extract_visium_assets:
    input:
        sp_input = sp_input
    output:
        img = f"{_assets_dir}/tissue_image.png",
        pos = f"{_assets_dir}/tissue_positions.csv",
        sf  = f"{_assets_dir}/scale_factor.txt"
    params:
        out_dir = _assets_dir
    singularity:
        _bayes_sif
    shell:
        "Rscript subworkflows/visualization/extract_visium_assets.R "
        "--sp_input {input.sp_input} --out_dir {params.out_dir}"

rule bayesspace_clustering:
    input:
        sp_input = sp_input
    output:
        _clustering
    params:
        q = str(_bayes_q)
    singularity:
        _bayes_sif
    shell:
        "Rscript subworkflows/clustering/bayesspace_cluster.R "
        "--sp_input {input.sp_input} --output {output} --q {params.q}"

rule seurat_clustering:
    input:
        sp_input = sp_input
    output:
        _seurat_cl
    params:
        res = str(_seurat_res)
    singularity:
        _bayes_sif
    shell:
        "Rscript subworkflows/clustering/seurat_cluster.R "
        "--sp_input {input.sp_input} --output {output} --resolution {params.res}"

rule generate_visualization:
    input:
        props      = _props,
        bayesspace = _clustering,
        seurat     = _seurat_cl,
        img        = rules.extract_visium_assets.output.img,
        pos        = rules.extract_visium_assets.output.pos,
        sf         = rules.extract_visium_assets.output.sf
    output:
        _html
    params:
        props_csv   = ",".join(_props),
        methods_csv = ",".join(methods),
        n_largest   = _n_largest
    singularity:
        _visu_sif
    shell:
        # Two clusterings -> the Seurat <-> BayesSpace toggle: primary = Seurat (transcriptomic),
        # secondary = BayesSpace (spatial).
        """
        SF=$(cat {input.sf})
        python3 subworkflows/visualization/sp_visualizer.py none \
            {params.props_csv} {input.pos} {input.seurat} {input.img} \
            {params.n_largest} "$SF" {output} {params.methods_csv} \
            {input.bayesspace} Seurat,BayesSpace
        """

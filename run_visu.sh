#!/usr/bin/env bash

################################ Slurm options #################################
#SBATCH --job-name=deconvolista_visu
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
# --output / --error / --partition: passed by launch.sh from config.env.
#
# Visualization-only job: builds the interactive HTML over ALL methods. The proportions are
# produced by the deconvolution jobs (CPU batch + cell2location/ddls on GPU); run this job AFTER
# them. run_all.sh submits it with an `--dependency=afterok:...` so it waits for the proportions,
# then Snakemake skips the (already done) deconvolution and only runs clustering + visualization.
################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'Job:' $SLURM_JOB_NAME '(' $SLURM_JOB_ID ')'
echo 'Node:' $HOSTNAME
echo '########################################'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CONFIG_FILE:-$SCRIPT_DIR/config.env}"

source "$CONDA_SETUP"
conda activate "$CONDA_ENV"

cd "$DECONVOLISTA_DIR"

CONFIG=(
    mode="run_dataset"
    methods="${VISU_METHODS:-rctd,nnls,spatialdwls,cell2location,ddls}"
    sc_input="${SC_INPUT:-unit-test/test_sc_data.rds}"
    sp_input="${SP_INPUT:-unit-test/test_sp_data.rds}"
    output="$OUTPUT_DIR"
    sif_dir="$SIF_DIR"
    use_gpu="false"
    skip_metrics="true"
    annot="${ANNOT:-subclass}"
    map_genes="false"
    do_visu="true"
    bayes_q="${BAYES_Q:-auto}"
)

snakemake -s main.smk --unlock --config "${CONFIG[@]}" || true

# --rerun-triggers mtime: rely on file timestamps only, so the already-computed proportions are
# NOT re-run (in particular cell2location is never accidentally re-run here on CPU).
snakemake -s main.smk -c 8 \
    --config "${CONFIG[@]}" \
    --rerun-triggers mtime \
    --use-singularity \
    --singularity-args "--bind $DATA_BIND --bind $OUTPUT_DIR" \
    --keep-going

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)

#!/usr/bin/env bash

################################ Slurm options #################################
#SBATCH --job-name=deconvolista
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=128G
# --output / --error / --partition: passed by launch.sh from config.env
################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'Job:' $SLURM_JOB_NAME '(' $SLURM_JOB_ID ')'
echo 'Node:' $HOSTNAME
echo '########################################'

# Load environment configuration (paths, conda). CONFIG_FILE is passed by launch.sh;
# otherwise fall back to config.env next to this script.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CONFIG_FILE:-$SCRIPT_DIR/config.env}"

source "$CONDA_SETUP"
conda activate "$CONDA_ENV"

cd "$DECONVOLISTA_DIR"

# Shared config, reused by the auto-unlock and the run
CONFIG=(
    mode="run_dataset"
    methods="${METHODS:-rctd,cell2location,nnls,spatialdwls,dirichlet,ddls}"
    sc_input="${SC_INPUT:-unit-test/test_sc_data.rds}"
    sp_input="${SP_INPUT:-unit-test/test_sp_data.rds}"
    output="$OUTPUT_DIR"
    sif_dir="$SIF_DIR"
    use_gpu="false"
    skip_metrics="false"
    annot="${ANNOT:-subclass}"
    map_genes="false"
    load_model="false"
)

# Auto-unlock (sequential jobs): avoids LockException after a scancel
snakemake -s main.smk --unlock --config "${CONFIG[@]}" || true

snakemake -s main.smk -c 12 \
    --config "${CONFIG[@]}" \
    --use-singularity \
    --singularity-args "--bind $DATA_BIND --bind $OUTPUT_DIR" \
    --keep-going

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)

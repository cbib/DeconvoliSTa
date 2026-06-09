#!/usr/bin/env bash

################################ Slurm options #################################
#SBATCH --job-name=c2l_gpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=100G
# --output / --error / --partition / --gres: passed by launch.sh from config.env.
# GPU_GRES defaults to a 1g.24gb H100 MIG (enough for the test data). For large
# real datasets, bump it in config.env (4g.47gb or full nvidia_h100_nvl).
################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'Job:' $SLURM_JOB_NAME '(' $SLURM_JOB_ID ')'
echo 'Node:' $HOSTNAME
echo '########################################'

# Load environment configuration (paths, conda, CUDA module)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CONFIG_FILE:-$SCRIPT_DIR/config.env}"

# CUDA host driver (so --nv can expose the GPU to the container)
[[ -n "${CUDA_MODULE:-}" ]] && module load "$CUDA_MODULE"

source "$CONDA_SETUP"
conda activate "$CONDA_ENV"

cd "$DECONVOLISTA_DIR"

# Sanity check: is the GPU visible on the host before starting Snakemake?
nvidia-smi || echo "WARNING: nvidia-smi failed - no GPU visible on host"

CONFIG=(
    mode="run_dataset"
    methods="cell2location"
    sc_input="${SC_INPUT:-unit-test/test_sc_data.rds}"
    sp_input="${SP_INPUT:-unit-test/test_sp_data.rds}"
    output="$OUTPUT_DIR"
    sif_dir="$SIF_DIR"
    use_gpu="true"
    skip_metrics="true"
    annot="${ANNOT:-subclass}"
    map_genes="false"
    load_model="false"
)

snakemake -s main.smk --unlock --config "${CONFIG[@]}" || true

snakemake -s main.smk -c 8 \
    --config "${CONFIG[@]}" \
    --use-singularity \
    --singularity-args "--nv --bind $DATA_BIND --bind $OUTPUT_DIR" \
    --keep-going

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)

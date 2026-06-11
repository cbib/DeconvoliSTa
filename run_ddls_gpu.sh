#!/usr/bin/env bash

################################ Slurm options #################################
#SBATCH --job-name=ddls_gpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
# --output / --error / --partition / --gres: passed by launch.sh from config.env.
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

# Sanity check: is the GPU visible on the host?
nvidia-smi || echo "WARNING: nvidia-smi failed - no GPU visible on host"

CONFIG=(
    mode="run_dataset"
    methods="ddls"
    sc_input="${SC_INPUT:-unit-test/test_sc_data.rds}"
    sp_input="${SP_INPUT:-unit-test/test_sp_data.rds}"
    output="$OUTPUT_DIR"
    sif_dir="$SIF_DIR"
    use_gpu="true"
    skip_metrics="true"
    annot="${ANNOT:-subclass}"
)

# Auto-unlock (sequential jobs)
snakemake -s main.smk --unlock --config "${CONFIG[@]}" || true

# NOTE: TensorFlow 2.8 ships no precompiled kernels for the H100 (compute 9.0)
# -> PTX JIT compilation on the first run (~30 min), then cached. Training then
# runs on GPU. If too slow/unstable, rebuild the DDLS image with TF >= 2.15.
snakemake -s main.smk -c 8 \
    --config "${CONFIG[@]}" \
    --use-singularity \
    --singularity-args "--nv --bind $DATA_BIND --bind $OUTPUT_DIR" \
    --keep-going

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)

#!/usr/bin/env bash

################################ Slurm options #################################

#SBATCH --job-name=c2l_gpu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:nvidia_h100_nvl_4g.47gb
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=100G
#SBATCH --output=/scratch/nmoualhi/c2l_gpu_%j.out
#SBATCH --error=/scratch/nmoualhi/c2l_gpu_%j.err

################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'Job:' $SLURM_JOB_NAME '(' $SLURM_JOB_ID ')'
echo 'Node:' $HOSTNAME
echo '########################################'

# CUDA host driver (needed so --nv can expose the GPU to the container).
# Verify the exact module name once with: module avail cuda
module load nvidia/cuda/12.9

source /isilon/modules/apps/conda/etc/profile.d/conda.sh
conda activate snakemake

cd /scratch/nmoualhi/DeconvoliSTa

# Sanity check: is the GPU visible on the host before we even start Snakemake?
nvidia-smi || echo "WARNING: nvidia-smi failed - no GPU visible on host"

snakemake -s main.smk -c 8 \
    --config \
    mode="run_dataset" \
    methods="cell2location" \
    sc_input="unit-test/test_sc_data.rds" \
    sp_input="unit-test/test_sp_data.rds" \
    output="/scratch/nmoualhi/res" \
    use_gpu="true" \
    skip_metrics="true" \
    annot="subclass" \
    map_genes="false" \
    load_model="false" \
    --use-singularity \
    --singularity-args '--nv --bind /mnt/cbib/RetinRNA/spatial' \
    --keep-going

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)

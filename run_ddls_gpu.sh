#!/usr/bin/env bash

################################ Slurm options #################################

#SBATCH --job-name=ddls_gpu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:nvidia_h100_nvl_1g.24gb
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --output=/scratch/nmoualhi/ddls_gpu_%j.out
#SBATCH --error=/scratch/nmoualhi/ddls_gpu_%j.err

################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'Job:' $SLURM_JOB_NAME '(' $SLURM_JOB_ID ')'
echo 'Node:' $HOSTNAME
echo '########################################'

# CUDA host driver (pour que --nv expose le GPU au container)
module load nvidia/cuda/12.9

source /isilon/modules/apps/conda/etc/profile.d/conda.sh
conda activate snakemake

cd /scratch/nmoualhi/DeconvoliSTa

# GPU visible sur l'hôte ?
nvidia-smi || echo "WARNING: nvidia-smi failed - no GPU visible on host"

CONFIG=(
    mode="run_dataset"
    methods="ddls"
    sc_input="unit-test/test_sc_data.rds"
    sp_input="unit-test/test_sp_data.rds"
    output="/scratch/nmoualhi/res"
    use_gpu="true"
    skip_metrics="true"
    annot="subclass"
)

# Déverrouillage auto (jobs séquentiels)
snakemake -s main.smk --unlock --config "${CONFIG[@]}" || true

# NB: TensorFlow 2.8 n'a pas de kernels precompiles pour la H100 (compute 9.0)
# -> JIT depuis PTX au premier run (~30 min), puis mis en cache. L'entrainement
# tourne ensuite sur GPU. Si trop lent/instable -> rebuild image avec TF >=2.15.
snakemake -s main.smk -c 8 \
    --config "${CONFIG[@]}" \
    --use-singularity \
    --singularity-args '--nv --bind /mnt/cbib/RetinRNA/spatial --bind /scratch/nmoualhi' \
    --keep-going

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)

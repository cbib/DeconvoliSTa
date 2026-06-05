#!/usr/bin/env bash

################################ Slurm options #################################

#SBATCH --job-name=deconvolista
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=128G
#SBATCH --output=/scratch/nmoualhi/deconvolista_%j.out
#SBATCH --error=/scratch/nmoualhi/deconvolista_%j.err

################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'Job:' $SLURM_JOB_NAME '(' $SLURM_JOB_ID ')'
echo 'Node:' $HOSTNAME
echo '########################################'

source /isilon/modules/apps/conda/etc/profile.d/conda.sh
conda activate snakemake

cd /scratch/nmoualhi/DeconvoliSTa

# Config commune, réutilisée par le déverrouillage auto et le run
CONFIG=(
    mode="run_dataset"
    methods="${METHODS:-rctd,cell2location,nnls,spatialdwls,dirichlet,ddls}"
    sc_input="unit-test/test_sc_data.rds"
    sp_input="unit-test/test_sp_data.rds"
    output="/scratch/nmoualhi/res"
    use_gpu="false"
    skip_metrics="false"
    annot="subclass"
    map_genes="false"
    load_model="false"
)

# Déverrouillage auto (jobs séquentiels, un seul à la fois) : évite le LockException après un scancel
snakemake -s main.smk --unlock --config "${CONFIG[@]}" || true

snakemake -s main.smk -c 12 \
    --config "${CONFIG[@]}" \
    --use-singularity \
    --singularity-args '--bind /mnt/cbib/RetinRNA/spatial --bind /scratch/nmoualhi' \
    --keep-going

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
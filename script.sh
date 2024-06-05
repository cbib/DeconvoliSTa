#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=demo_job

### Limit run time "hours:minutes:seconds" (default: 365 days)
#SBATCH --time=01:00:00

### Specify requirements - Task (default: 1 node, 1 Core, 12.5G mem/cpu)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=12800MB

################################################################################

# Useful information to print
echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
# Detail Information:
scontrol show job $SLURM_JOB_ID
echo '########################################'

# modules loading
module add python


# What you aco lauh. Awesome.'

nextflow run main.nf -profile local,docker --methods music --sc_input unit-test/test_sc_data.rds \
--sp_input unit-test/test_sp_data.rds --annot subclass

echo '###########'

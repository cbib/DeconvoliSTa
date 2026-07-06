#!/usr/bin/env bash
#SBATCH --job-name=deconvolista
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=64G
# NB: --job-name / --time / --mem / --partition / --gres / --output / --error are all set
# per-submission by ./deconvolista (sbatch CLI flags override these #SBATCH defaults). The
# values above are only fallbacks for a manual `sbatch slurm_job.sh`.

# =============================================================================
# Generic DeconvoliSTa SLURM job body — ONE script for every case (CPU methods, a GPU
# method, or the visualization). The variation is driven entirely by these env vars,
# exported by ./deconvolista:
#   RUN_METHODS   comma-separated methods to run (e.g. "nnls,rctd" or "cell2location")
#   RUN_USE_GPU   "true"  -> expose the GPU to the container (--nv) + load the CUDA module
#   RUN_DO_VISU   "true"  -> also build the interactive visualization (clustering + HTML)
# Data inputs (SC_INPUT / SP_INPUT / ANNOT / BAYES_Q ...) come from config.env, overridable.
# =============================================================================
set -uo pipefail

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'Job:' "${SLURM_JOB_NAME:-local}" '(' "${SLURM_JOB_ID:-?}" ')'
echo 'Node:' "$HOSTNAME"
echo 'Methods:' "${RUN_METHODS:?RUN_METHODS not set — launch via ./deconvolista}" \
     '| GPU:' "${RUN_USE_GPU:-false}" '| visu:' "${RUN_DO_VISU:-false}"
echo '########################################'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CONFIG_FILE:-$SCRIPT_DIR/config.env}"

source "$CONDA_SETUP"
conda activate "$CONDA_ENV"
cd "$DECONVOLISTA_DIR"

# Redirect Singularity's cache + build tmpdir onto the big filesystem (OUTPUT_DIR lives there):
# pulling/converting the large GPU images (cell2location ~9 GB, ddls ~10 GB) to SIF overflows the
# default ~/.singularity cache and /tmp on HPC compute nodes ("No space left on device").
export APPTAINER_CACHEDIR="${SINGULARITY_CACHEDIR:-$OUTPUT_DIR/.singularity/cache}"
export APPTAINER_TMPDIR="${SINGULARITY_TMPDIR:-$OUTPUT_DIR/.singularity/tmp}"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

# Singularity binds shared by every job.
sing_args="--bind $DATA_BIND --bind $OUTPUT_DIR"

# GPU jobs: load the host CUDA driver and expose the GPU to the container (--nv).
if [[ "${RUN_USE_GPU:-false}" == "true" ]]; then
    [[ -n "${CUDA_MODULE:-}" ]] && module load "$CUDA_MODULE"
    nvidia-smi || echo "WARNING: nvidia-smi failed — no GPU visible on host"
    sing_args="--nv $sing_args"
fi

CONFIG=(
    mode="run_dataset"
    methods="$RUN_METHODS"
    sc_input="${SC_INPUT:-unit-test/test_sc_data.rds}"
    sp_input="${SP_INPUT:-unit-test/test_sp_data.rds}"
    output="$OUTPUT_DIR"
    sif_dir="$SIF_DIR"
    use_gpu="${RUN_USE_GPU:-false}"
    skip_metrics="${SKIP_METRICS:-true}"
    annot="${ANNOT:-subclass}"
    map_genes="${MAP_GENES:-false}"
    load_model="false"
    do_visu="${RUN_DO_VISU:-false}"
    bayes_q="${BAYES_Q:-auto}"
)

# Auto-unlock (in case a previous job was cancelled mid-run, leaving a stale lock).
snakemake -s main.smk --unlock --config "${CONFIG[@]}" || true

# The visualization job reuses proportions already produced by the deconvolution jobs;
# --rerun-triggers mtime keeps Snakemake from re-running them (esp. cell2location on CPU).
rerun_args=()
[[ "${RUN_DO_VISU:-false}" == "true" ]] && rerun_args=(--rerun-triggers mtime)

snakemake -s main.smk -c "${SLURM_CPUS_ON_NODE:-8}" \
    --config "${CONFIG[@]}" \
    --use-singularity \
    --singularity-args "$sing_args" \
    "${rerun_args[@]}" \
    --keep-going
status=$?   # propagate Snakemake's exit code so SLURM marks the job failed (the visu's
            # afterok dependency then won't run on incomplete data).

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
exit $status

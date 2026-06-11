#!/usr/bin/env bash
# =============================================================================
# DeconvoliSTa launch wrapper: reads config.env and submits the right SLURM job
# (partition, GPU, log directories) with no hardcoded paths in the scripts.
#
# Usage:
#   ./launch.sh run_deconvolista.sh                          # CPU batch (all methods)
#   ./launch.sh run_cell2location_gpu.sh                     # cell2location on GPU
#   ./launch.sh run_ddls_gpu.sh                              # DDLS on GPU
#   METHODS=nnls,dirichlet ./launch.sh run_deconvolista.sh   # on-the-fly override
#
# Prerequisite: cp config.env.example config.env && nano config.env
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CONFIG_FILE="${CONFIG_FILE:-$SCRIPT_DIR/config.env}"

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "ERROR: $CONFIG_FILE not found." >&2
  echo "  -> cp config.env.example config.env  then edit it for your environment." >&2
  exit 1
fi
source "$CONFIG_FILE"

script="${1:?usage: ./launch.sh <run_*.sh>  (e.g. run_deconvolista.sh)}"
[[ -f "$SCRIPT_DIR/$script" ]] || { echo "ERROR: script not found: $script" >&2; exit 1; }

mkdir -p "$LOG_DIR" 2>/dev/null || true

# --export=ALL forwards the current environment (CONFIG_FILE, METHODS, SC_INPUT...) to the job.
sbatch_opts=(--output="$LOG_DIR/%x_%j.out" --error="$LOG_DIR/%x_%j.err" --export=ALL)

if [[ "$script" == *gpu* ]]; then
  echo "GPU submission: $script  (partition=$PARTITION_GPU, gres=$GPU_GRES)"
  sbatch "${sbatch_opts[@]}" --partition="$PARTITION_GPU" --gres="$GPU_GRES" "$SCRIPT_DIR/$script"
else
  echo "CPU submission: $script  (partition=$PARTITION_CPU)"
  sbatch "${sbatch_opts[@]}" --partition="$PARTITION_CPU" "$SCRIPT_DIR/$script"
fi

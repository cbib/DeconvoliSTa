#!/usr/bin/env bash
# =============================================================================
# DeconvoliSTa launch wrapper: reads config.env and submits the right SLURM job
# (partition, GPU, log directories) with no hardcoded paths in the scripts.
#
# Usage:
#   ./launch.sh run_deconvolista.sh                          # CPU batch (CPU methods)
#   ./launch.sh run_cell2location_gpu.sh                     # cell2location on GPU
#   ./launch.sh run_ddls_gpu.sh                              # DDLS on GPU
#   METHODS=nnls,dirichlet ./launch.sh run_deconvolista.sh   # on-the-fly override
#   ./launch.sh --dependency=afterok:123 run_visu.sh         # extra sbatch flags (passed through)
#
# Prints the submitted job ID on stdout (so jobs can be chained); progress messages go to stderr.
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

# Optional sbatch flags (e.g. --dependency=afterok:123) given before the script name.
passthrough=()
while [[ "${1:-}" == --* ]]; do passthrough+=("$1"); shift; done

script="${1:?usage: ./launch.sh [--sbatch-flags] <run_*.sh>  (e.g. run_deconvolista.sh)}"
[[ -f "$SCRIPT_DIR/$script" ]] || { echo "ERROR: script not found: $script" >&2; exit 1; }

mkdir -p "$LOG_DIR" 2>/dev/null || true

# --parsable -> sbatch prints only the job ID (captured by run_all.sh to chain dependencies).
# --export=ALL forwards the current environment (CONFIG_FILE, METHODS, SC_INPUT...) to the job.
sbatch_opts=(--parsable --output="$LOG_DIR/%x_%j.out" --error="$LOG_DIR/%x_%j.err" --export=ALL "${passthrough[@]}")

if [[ "$script" == *gpu* ]]; then
  echo "GPU submission: $script  (partition=$PARTITION_GPU, gres=$GPU_GRES)" >&2
  sbatch "${sbatch_opts[@]}" --partition="$PARTITION_GPU" --gres="$GPU_GRES" "$SCRIPT_DIR/$script"
else
  echo "CPU submission: $script  (partition=$PARTITION_CPU)" >&2
  sbatch "${sbatch_opts[@]}" --partition="$PARTITION_CPU" "$SCRIPT_DIR/$script"
fi

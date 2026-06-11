#!/usr/bin/env bash
# =============================================================================
# Run the FULL pipeline in one command: the CPU methods on the CPU partition and
# the GPU methods (cell2location, ddls) on the GPU partition -- each on the right
# hardware, as three SLURM jobs submitted at once.
#
# Usage:
#   ./run_all.sh
#   SC_INPUT=ref.rds SP_INPUT=sample.rds ANNOT=cell_type ./run_all.sh
#
# Data/annotation overrides (SC_INPUT, SP_INPUT, ANNOT, ...) are forwarded to every
# job by launch.sh (--export=ALL). Prerequisite: cp config.env.example config.env.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"$SCRIPT_DIR/launch.sh" run_deconvolista.sh        # rctd, nnls, spatialdwls, dirichlet (CPU)
"$SCRIPT_DIR/launch.sh" run_cell2location_gpu.sh   # cell2location (GPU)
"$SCRIPT_DIR/launch.sh" run_ddls_gpu.sh            # ddls (GPU)

echo "Submitted: CPU batch + cell2location (GPU) + ddls (GPU). Track with: squeue -u \$USER"

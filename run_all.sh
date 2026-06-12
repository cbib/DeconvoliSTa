#!/usr/bin/env bash
# =============================================================================
# Run the pipeline on SLURM, your way: pick any set of methods and run_all routes each one to the
# right hardware (CPU vs GPU), then (optionally) builds the interactive visualization over exactly
# the methods you ran — with a SLURM dependency so the visu waits for the deconvolution.
#
# Usage:
#   SC_INPUT=ref.rds SP_INPUT=sample.rds ANNOT=cell_type ./run_all.sh
#   METHODS=cell2location ...            ./run_all.sh    # just one method
#   METHODS=rctd,cell2location ...       ./run_all.sh    # any subset
#   DO_VISU=false METHODS=rctd ...       ./run_all.sh    # proportions only, no visualization
#
# Env vars (all optional):
#   METHODS   default "rctd,nnls,spatialdwls,cell2location,ddls"  (any subset; routed automatically)
#   DO_VISU   default "true"   ("false" = deconvolution only)
#   BAYES_Q   BayesSpace q for the visualization clustering (default "auto")
#
# Routing: rctd / nnls / spatialdwls / dirichlet -> CPU partition (one job);
#          cell2location, ddls -> GPU partition (one job each).
#
# Prerequisite: cp config.env.example config.env. For DO_VISU=true, build sp_bayesspace.sif and
# visu.sif once (see SETUP.md).
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

METHODS="${METHODS:-rctd,nnls,spatialdwls,cell2location,ddls}"

# --- classify the requested methods into CPU vs GPU ---
cpu_methods=(); want_c2l=0; want_ddls=0
IFS=',' read -ra _arr <<< "$METHODS"
for m in "${_arr[@]}"; do
  case "$m" in
    cell2location) want_c2l=1 ;;
    ddls)          want_ddls=1 ;;
    rctd|nnls|spatialdwls|dirichlet) cpu_methods+=("$m") ;;
    *) echo "ERROR: unknown method '$m' (use rctd,nnls,spatialdwls,dirichlet,cell2location,ddls)" >&2; exit 1 ;;
  esac
done

deps=()

# --- CPU job (the CPU methods, if any) — no visualization here ---
if [[ ${#cpu_methods[@]} -gt 0 ]]; then
  cpu_csv=$(IFS=,; echo "${cpu_methods[*]}")
  cpu_id=$(DO_VISU=false METHODS="$cpu_csv" "$SCRIPT_DIR/launch.sh" run_deconvolista.sh)
  deps+=("$cpu_id")
  echo "  CPU ($cpu_csv): job $cpu_id" >&2
fi

# --- GPU jobs (one per GPU method) ---
if [[ $want_c2l -eq 1 ]]; then
  c2l_id=$("$SCRIPT_DIR/launch.sh" run_cell2location_gpu.sh); deps+=("$c2l_id")
  echo "  cell2location (GPU): job $c2l_id" >&2
fi
if [[ $want_ddls -eq 1 ]]; then
  ddls_id=$("$SCRIPT_DIR/launch.sh" run_ddls_gpu.sh); deps+=("$ddls_id")
  echo "  ddls (GPU): job $ddls_id" >&2
fi

# --- Visualization job, waiting for all the deconvolution jobs (unless DO_VISU=false) ---
if [[ "${DO_VISU:-true}" == "true" ]]; then
  dep=$(IFS=:; echo "afterok:${deps[*]}")
  visu_id=$(VISU_METHODS="$METHODS" "$SCRIPT_DIR/launch.sh" --dependency="$dep" run_visu.sh)
  echo "  visualization (after the above): job $visu_id" >&2
fi

echo "Submitted methods: $METHODS  (visualization: ${DO_VISU:-true}). Track with: squeue -u \$USER"

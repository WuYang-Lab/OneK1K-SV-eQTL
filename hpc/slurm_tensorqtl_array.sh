#!/usr/bin/env bash
#SBATCH --job-name=sv_eqtl
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=72:00:00
#SBATCH --gres=gpu:1
#SBATCH --output=logs/slurm/tensorqtl_%A_%a.out
#SBATCH --error=logs/slurm/tensorqtl_%A_%a.err

set -Eeuo pipefail

CONFIG="${1:?Usage: sbatch --array=1-N hpc/slurm_tensorqtl_array.sh <config.env> [sv|vntr|joint]}"
GENOTYPE_LABEL="${2:-sv}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# shellcheck disable=SC1090
source "$CONFIG"
: "${CELL_TYPES_FILE:?CELL_TYPES_FILE is not set in the configuration}"
mapfile -t CELL_TYPES < <(awk 'NF && $1 !~ /^#/ {print $1}' "$CELL_TYPES_FILE")
INDEX=$((SLURM_ARRAY_TASK_ID - 1))
(( INDEX >= 0 && INDEX < ${#CELL_TYPES[@]} )) || {
  printf '[ERROR] SLURM_ARRAY_TASK_ID is outside the cell-type list\n' >&2
  exit 2
}
CELL_TYPE="${CELL_TYPES[$INDEX]}"

bash "${REPO_ROOT}/scripts/06_run_tensorqtl.sh" "$CONFIG" "$CELL_TYPE" "$GENOTYPE_LABEL"

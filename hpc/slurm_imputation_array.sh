#!/usr/bin/env bash
#SBATCH --job-name=sv_impute
#SBATCH --array=1-22
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=72:00:00
#SBATCH --output=logs/slurm/impute_%A_%a.out
#SBATCH --error=logs/slurm/impute_%A_%a.err

set -Eeuo pipefail

CONFIG="${1:?Usage: sbatch hpc/slurm_imputation_array.sh <config.env> [SV|VNTR]}"
PANEL="${2:-SV}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mapfile -t CHROMOSOMES < <(awk 'NF && $1 !~ /^#/ {print $1}' "${REPO_ROOT}/config/chromosomes.txt")
CHR="${CHROMOSOMES[$((SLURM_ARRAY_TASK_ID - 1))]}"

bash "${REPO_ROOT}/scripts/03_impute_chromosome.sh" "$CONFIG" "$PANEL" "$CHR"

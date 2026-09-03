#!/usr/bin/env bash
#SBATCH --job-name=sv_qc
#SBATCH --array=1-22
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/qc_%A_%a.out
#SBATCH --error=logs/slurm/qc_%A_%a.err

set -Eeuo pipefail

CONFIG="${1:?Usage: sbatch hpc/slurm_filter_array.sh <config.env> [SV|VNTR]}"
PANEL="${2:-SV}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mapfile -t CHROMOSOMES < <(awk 'NF && $1 !~ /^#/ {print $1}' "${REPO_ROOT}/config/chromosomes.txt")
CHR="${CHROMOSOMES[$((SLURM_ARRAY_TASK_ID - 1))]}"

bash "${REPO_ROOT}/scripts/04_filter_imputed_chromosome.sh" "$CONFIG" "$PANEL" "$CHR"

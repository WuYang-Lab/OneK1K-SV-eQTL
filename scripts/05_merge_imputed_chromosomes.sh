#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

load_config "${1:-}"
PANEL="${2:-}"
PANEL_UPPER=$(printf '%s' "$PANEL" | tr '[:lower:]' '[:upper:]')
[[ "$PANEL_UPPER" == "SV" || "$PANEL_UPPER" == "VNTR" ]] || die "Panel must be SV or VNTR"
PANEL_LOWER=$(printf '%s' "$PANEL_UPPER" | tr '[:upper:]' '[:lower:]')
require_var CHROMOSOMES_FILE
require_command plink2
require_file "$CHROMOSOMES_FILE"

OUT_DIR="${PROJECT_ROOT}/results/genotypes"
LOG_DIR="${PROJECT_ROOT}/logs/qc/${PANEL_LOWER}"
MERGE_LIST="${OUT_DIR}/${PANEL_LOWER}.pmerge_list.txt"
PGEN_PREFIX="${OUT_DIR}/${PANEL_LOWER}.merged"
BED_PREFIX="${OUT_DIR}/${PANEL_LOWER}.merged.bedfmt"
mkdir -p "$OUT_DIR" "$LOG_DIR"
: > "$MERGE_LIST"

while read -r CHR_NUM; do
  PREFIX="${PROJECT_ROOT}/results/qc/${PANEL_LOWER}/chr${CHR_NUM}/target.hg38.chr${CHR_NUM}.${PANEL_UPPER}_imputed.qc"
  require_file "${PREFIX}.pgen"
  require_file "${PREFIX}.pvar"
  require_file "${PREFIX}.psam"
  printf '%s\n' "$PREFIX" >> "$MERGE_LIST"
done < <(clean_list "$CHROMOSOMES_FILE")

log "Merging chromosome-wise ${PANEL_UPPER} PGEN files"
plink2 \
  --pmerge-list "$MERGE_LIST" \
  --make-pgen \
  --out "$PGEN_PREFIX" \
  > "${LOG_DIR}/merge_pgen.log" 2>&1

plink2 \
  --pfile "$PGEN_PREFIX" \
  --make-bed \
  --out "$BED_PREFIX" \
  > "${LOG_DIR}/merge_make_bed.log" 2>&1

VARIANT_N=$(awk 'END {print NR-1}' "${PGEN_PREFIX}.pvar")
SAMPLE_N=$(awk 'END {print NR-1}' "${PGEN_PREFIX}.psam")
{
  printf 'panel\tvariants\tsamples\n'
  printf '%s\t%s\t%s\n' "$PANEL_UPPER" "$VARIANT_N" "$SAMPLE_N"
} > "${OUT_DIR}/${PANEL_LOWER}.merge_summary.tsv"

log "Merged ${PANEL_UPPER} set: ${VARIANT_N} variants and ${SAMPLE_N} samples"
log "PGEN prefix: ${PGEN_PREFIX}"
log "BED prefix: ${BED_PREFIX}"

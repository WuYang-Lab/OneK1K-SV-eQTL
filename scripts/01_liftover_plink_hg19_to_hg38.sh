#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

load_config "${1:-}"
require_var BFILE_HG19
require_var LIFTOVER_BIN
require_var HG19_TO_HG38_CHAIN
require_var CHROMOSOMES_FILE
require_command plink
require_command awk
require_command sort
require_file "${BFILE_HG19}.bed"
require_file "${BFILE_HG19}.bim"
require_file "${BFILE_HG19}.fam"
require_file "$HG19_TO_HG38_CHAIN"
require_file "$CHROMOSOMES_FILE"

if [[ "$LIFTOVER_BIN" == */* ]]; then
  require_file "$LIFTOVER_BIN"
else
  require_command "$LIFTOVER_BIN"
fi

OUT_DIR="${PROJECT_ROOT}/results/genotypes"
WORK_DIR="${PROJECT_ROOT}/work/liftover"
LOG_DIR="${PROJECT_ROOT}/logs"
OUT_PREFIX="${OUT_DIR}/target.hg38"
REMOVE_DUPLICATE_POSITIONS="${REMOVE_DUPLICATE_POSITIONS:-true}"
mkdir -p "$OUT_DIR" "$WORK_DIR" "$LOG_DIR"

INPUT_ID_DUP="${WORK_DIR}/input.duplicate_variant_ids.txt"
awk '{count[$2]++} END {for (id in count) if (count[id] > 1) print id}' \
  "${BFILE_HG19}.bim" | sort > "$INPUT_ID_DUP"
[[ ! -s "$INPUT_ID_DUP" ]] || die "Input BIM has duplicate variant IDs; see ${INPUT_ID_DUP}"

HG19_BED="${WORK_DIR}/target.hg19.bed"
clean_list "$CHROMOSOMES_FILE" > "${WORK_DIR}/chromosomes.clean.txt"
awk 'BEGIN{OFS="\t"}
  NR==FNR {gsub(/^chr/, "", $1); keep[$1]=1; next}
  {chrom=$1; gsub(/^chr/, "", chrom)}
  chrom in keep {print "chr" chrom, $4-1, $4, $2}' \
  "${WORK_DIR}/chromosomes.clean.txt" "${BFILE_HG19}.bim" > "$HG19_BED"
require_file "$HG19_BED"

log "Lifting PLINK variant coordinates from hg19 to hg38"
"$LIFTOVER_BIN" \
  "$HG19_BED" \
  "$HG19_TO_HG38_CHAIN" \
  "${WORK_DIR}/target.hg38.bed" \
  "${WORK_DIR}/target.hg19_to_hg38.unmapped.bed"

awk 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {
  chrom=$1; sub(/^chr/, "", chrom); print $4, chrom
}' "${WORK_DIR}/target.hg38.bed" > "${WORK_DIR}/update_chr.tsv"

awk 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {print $4, $3}' \
  "${WORK_DIR}/target.hg38.bed" > "${WORK_DIR}/update_pos.tsv"

cut -f1 "${WORK_DIR}/update_chr.tsv" > "${WORK_DIR}/lifted_variant_ids.txt"
require_file "${WORK_DIR}/lifted_variant_ids.txt"

plink \
  --bfile "$BFILE_HG19" \
  --extract "${WORK_DIR}/lifted_variant_ids.txt" \
  --make-bed \
  --out "${WORK_DIR}/target.hg19.lifted_only" \
  > "${LOG_DIR}/01_plink_extract.log" 2>&1

plink \
  --bfile "${WORK_DIR}/target.hg19.lifted_only" \
  --update-chr "${WORK_DIR}/update_chr.tsv" 2 1 \
  --update-map "${WORK_DIR}/update_pos.tsv" 2 1 \
  --make-bed \
  --out "${WORK_DIR}/target.hg38.with_duplicates" \
  > "${LOG_DIR}/01_plink_update_coordinates.log" 2>&1

DUP_POS="${WORK_DIR}/hg38.duplicate_position_variant_ids.txt"
awk '{key=$1 FS $4; count[key]++; ids[key]=ids[key] $2 FS}
  END {for (key in count) if (count[key] > 1) {
    n=split(ids[key], x, FS); for (i=1; i<=n; i++) if (x[i] != "") print x[i]
  }}' "${WORK_DIR}/target.hg38.with_duplicates.bim" | sort -u > "$DUP_POS"

if [[ -s "$DUP_POS" && "$REMOVE_DUPLICATE_POSITIONS" == "true" ]]; then
  log "Removing variants at duplicated hg38 coordinates"
  plink \
    --bfile "${WORK_DIR}/target.hg38.with_duplicates" \
    --exclude "$DUP_POS" \
    --make-bed \
    --out "$OUT_PREFIX" \
    > "${LOG_DIR}/01_plink_remove_duplicate_positions.log" 2>&1
else
  plink \
    --bfile "${WORK_DIR}/target.hg38.with_duplicates" \
    --make-bed \
    --out "$OUT_PREFIX" \
    > "${LOG_DIR}/01_plink_finalize.log" 2>&1
fi

TOTAL=$(wc -l < "${BFILE_HG19}.bim")
LIFTED=$(wc -l < "${WORK_DIR}/lifted_variant_ids.txt")
FINAL=$(wc -l < "${OUT_PREFIX}.bim")
UNMAPPED=$(awk 'NF && $1 !~ /^#/ {n++} END {print n+0}' "${WORK_DIR}/target.hg19_to_hg38.unmapped.bed")
DUPLICATE_COORDINATES=$(wc -l < "$DUP_POS")

{
  printf 'metric\tvalue\n'
  printf 'input_variants\t%s\n' "$TOTAL"
  printf 'lifted_variants\t%s\n' "$LIFTED"
  printf 'unmapped_variants\t%s\n' "$UNMAPPED"
  printf 'variant_ids_at_duplicate_hg38_positions\t%s\n' "$DUPLICATE_COORDINATES"
  printf 'final_variants\t%s\n' "$FINAL"
} > "${OUT_DIR}/liftover_summary.tsv"

log "Completed liftOver: ${FINAL}/${TOTAL} variants retained"
log "Output prefix: ${OUT_PREFIX}"

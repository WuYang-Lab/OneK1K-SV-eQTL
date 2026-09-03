#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

load_config "${1:-}"
require_var SNP_BFILE_HG38
require_command plink

SV_BFILE="${SV_BFILE:-${PROJECT_ROOT}/results/genotypes/sv.merged.bedfmt}"
OUT_DIR="${PROJECT_ROOT}/results/genotypes/joint"
WORK_DIR="${PROJECT_ROOT}/work/joint_genotypes"
LOG_DIR="${PROJECT_ROOT}/logs"
OUT_PREFIX="${OUT_DIR}/snp_sv.joint"
mkdir -p "$OUT_DIR" "$WORK_DIR" "$LOG_DIR"

for suffix in bed bim fam; do
  require_file "${SNP_BFILE_HG38}.${suffix}"
  require_file "${SV_BFILE}.${suffix}"
done

awk '{print $1 "\t" $2}' "${SNP_BFILE_HG38}.fam" > "${WORK_DIR}/snp.samples.tsv"
awk '{print $1 "\t" $2}' "${SV_BFILE}.fam" > "${WORK_DIR}/sv.samples.tsv"
cmp -s "${WORK_DIR}/snp.samples.tsv" "${WORK_DIR}/sv.samples.tsv" \
  || die "SNP and SV FID/IID lists or their order differ"

awk '{print $2}' "${SNP_BFILE_HG38}.bim" | sort > "${WORK_DIR}/snp.variant_ids.txt"
awk '{print $2}' "${SV_BFILE}.bim" | sort > "${WORK_DIR}/sv.variant_ids.txt"
comm -12 "${WORK_DIR}/snp.variant_ids.txt" "${WORK_DIR}/sv.variant_ids.txt" \
  > "${WORK_DIR}/overlapping_variant_ids.txt"
[[ ! -s "${WORK_DIR}/overlapping_variant_ids.txt" ]] \
  || die "SNP and SV sets contain overlapping variant IDs; see ${WORK_DIR}/overlapping_variant_ids.txt"

log "Merging hg38 SNP and SV genotype sets"
if ! plink \
  --bfile "$SNP_BFILE_HG38" \
  --bmerge "${SV_BFILE}.bed" "${SV_BFILE}.bim" "${SV_BFILE}.fam" \
  --make-bed \
  --out "$OUT_PREFIX" \
  > "${LOG_DIR}/07_plink_merge_snp_sv.log" 2>&1; then
  die "PLINK SNP+SV merge failed; inspect ${LOG_DIR}/07_plink_merge_snp_sv.log and any .missnp file"
fi

SNP_N=$(wc -l < "${SNP_BFILE_HG38}.bim")
SV_N=$(wc -l < "${SV_BFILE}.bim")
JOINT_N=$(wc -l < "${OUT_PREFIX}.bim")
SAMPLE_N=$(wc -l < "${OUT_PREFIX}.fam")
{
  printf 'metric\tvalue\n'
  printf 'snp_variants\t%s\n' "$SNP_N"
  printf 'sv_variants\t%s\n' "$SV_N"
  printf 'joint_variants\t%s\n' "$JOINT_N"
  printf 'samples\t%s\n' "$SAMPLE_N"
} > "${OUT_DIR}/merge_summary.tsv"

log "Created joint genotype set with ${JOINT_N} variants and ${SAMPLE_N} samples"
log "Output prefix: ${OUT_PREFIX}"

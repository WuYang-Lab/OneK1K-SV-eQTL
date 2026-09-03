#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

load_config "${1:-}"
PANEL="${2:-}"
CHR_NUM=$(normalize_chr_number "${3:-}")
PANEL_UPPER=$(printf '%s' "$PANEL" | tr '[:lower:]' '[:upper:]')
[[ "$PANEL_UPPER" == "SV" || "$PANEL_UPPER" == "VNTR" ]] || die "Panel must be SV or VNTR"

PANEL_VAR="${PANEL_UPPER}_REFERENCE_VCF"
require_var "$PANEL_VAR"
REFERENCE_VCF="${!PANEL_VAR}"
require_var BEAGLE_JAR
require_command bcftools
require_command tabix
require_command java
require_command grep
require_file "$REFERENCE_VCF"
require_file "$BEAGLE_JAR"

TARGET_VCF="${TARGET_VCF:-${PROJECT_ROOT}/results/target/target.hg38.chr.id.vcf.gz}"
require_file "$TARGET_VCF"
CHR="chr${CHR_NUM}"
PANEL_LOWER=$(printf '%s' "$PANEL_UPPER" | tr '[:upper:]' '[:lower:]')
OUT_DIR="${PROJECT_ROOT}/results/imputation/${PANEL_LOWER}/${CHR}"
LOG_DIR="${PROJECT_ROOT}/logs/imputation/${PANEL_LOWER}"
mkdir -p "$OUT_DIR" "$LOG_DIR"

TARGET_CHR="${OUT_DIR}/target.hg38.${CHR}.vcf.gz"
REFERENCE_CHR="${OUT_DIR}/${PANEL_UPPER}_panel.hg38.${CHR}.vcf.gz"
OUT_PREFIX="${OUT_DIR}/target.hg38.${CHR}.${PANEL_UPPER}_imputed"
NTHREADS="${NTHREADS:-16}"
BEAGLE_MEMORY_GB="${BEAGLE_MEMORY_GB:-32}"

bcftools index -s "$TARGET_VCF" | cut -f1 | grep -Fxq "$CHR" \
  || die "Target VCF does not contain contig ${CHR}"
bcftools index -s "$REFERENCE_VCF" | cut -f1 | grep -Fxq "$CHR" \
  || die "Reference VCF does not contain contig ${CHR}"

log "Extracting target and ${PANEL_UPPER} reference records for ${CHR}"
bcftools view --regions "$CHR" -Oz -o "$TARGET_CHR" "$TARGET_VCF"
bcftools view --regions "$CHR" -Oz -o "$REFERENCE_CHR" "$REFERENCE_VCF"
tabix -f -p vcf "$TARGET_CHR"
tabix -f -p vcf "$REFERENCE_CHR"

TARGET_N=$(vcf_count "$TARGET_CHR")
REFERENCE_N=$(vcf_count "$REFERENCE_CHR")
(( TARGET_N > 0 )) || die "No target variants found on ${CHR}"
(( REFERENCE_N > 0 )) || die "No reference variants found on ${CHR}"

OVERLAP_N=$(bcftools isec -n=2 -w1 "$TARGET_CHR" "$REFERENCE_CHR" | awk '!/^#/ {n++} END {print n+0}')
log "${CHR}: target=${TARGET_N}, reference=${REFERENCE_N}, exact overlap=${OVERLAP_N}"

log "Running Beagle ${PANEL_UPPER} imputation for ${CHR}"
java "-Xmx${BEAGLE_MEMORY_GB}g" -jar "$BEAGLE_JAR" \
  "gt=${TARGET_CHR}" \
  "ref=${REFERENCE_CHR}" \
  "out=${OUT_PREFIX}" \
  "chrom=${CHR}" \
  "nthreads=${NTHREADS}" \
  > "${LOG_DIR}/${CHR}.stdout.log" 2>&1

IMPUTED_VCF="${OUT_PREFIX}.vcf.gz"
require_file "$IMPUTED_VCF"
tabix -f -p vcf "$IMPUTED_VCF"

HEADER=$(bcftools view -h "$IMPUTED_VCF")
grep -q 'ID=DR2' <<< "$HEADER" || die "Beagle output is missing INFO/DR2"
grep -q 'ID=AF' <<< "$HEADER" || die "Beagle output is missing INFO/AF"
grep -q 'ID=DS' <<< "$HEADER" || die "Beagle output is missing FORMAT/DS"
IMPUTED_N=$(vcf_count "$IMPUTED_VCF")

{
  printf 'chrom\tpanel\ttarget_variants\treference_variants\texact_overlap\timputed_variants\n'
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$CHR" "$PANEL_UPPER" "$TARGET_N" "$REFERENCE_N" "$OVERLAP_N" "$IMPUTED_N"
} > "${OUT_DIR}/imputation_counts.tsv"

log "Completed ${PANEL_UPPER} imputation for ${CHR}: ${IMPUTED_N} output variants"

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
PANEL_LOWER=$(printf '%s' "$PANEL_UPPER" | tr '[:upper:]' '[:lower:]')
CHR="chr${CHR_NUM}"

require_command bcftools
require_command tabix
require_command plink2

DR2_MIN="${DR2_MIN:-0.3}"
MAF_MIN="${MAF_MIN:-0.01}"
MAF_MAX="${MAF_MAX:-0.99}"
GENO_MAX="${GENO_MAX:-0.05}"
HWE_MIN="${HWE_MIN:-1e-6}"

IN_VCF="${PROJECT_ROOT}/results/imputation/${PANEL_LOWER}/${CHR}/target.hg38.${CHR}.${PANEL_UPPER}_imputed.vcf.gz"
OUT_DIR="${PROJECT_ROOT}/results/qc/${PANEL_LOWER}/${CHR}"
LOG_DIR="${PROJECT_ROOT}/logs/qc/${PANEL_LOWER}"
FILTERED_VCF="${OUT_DIR}/target.hg38.${CHR}.${PANEL_UPPER}_imputed.qc.vcf.gz"
PGEN_PREFIX="${OUT_DIR}/target.hg38.${CHR}.${PANEL_UPPER}_imputed.qc"
BED_PREFIX="${PGEN_PREFIX}.bedfmt"
mkdir -p "$OUT_DIR" "$LOG_DIR"
require_file "$IN_VCF"

BEFORE_N=$(vcf_count "$IN_VCF")
log "Filtering ${PANEL_UPPER} variants on ${CHR} with DR2>=${DR2_MIN} and AF in [${MAF_MIN}, ${MAF_MAX}]"
bcftools view \
  --include "INFO/DR2>=${DR2_MIN} && INFO/AF>=${MAF_MIN} && INFO/AF<=${MAF_MAX}" \
  --min-alleles 2 \
  --max-alleles 2 \
  -Oz \
  -o "$FILTERED_VCF" \
  "$IN_VCF"
tabix -f -p vcf "$FILTERED_VCF"
AF_DR2_N=$(vcf_count "$FILTERED_VCF")

plink2 \
  --vcf "$FILTERED_VCF" dosage=DS \
  --maf "$MAF_MIN" \
  --geno "$GENO_MAX" \
  --hwe "$HWE_MIN" \
  --make-pgen \
  --out "$PGEN_PREFIX" \
  > "${LOG_DIR}/${CHR}.plink2_qc.log" 2>&1

plink2 \
  --pfile "$PGEN_PREFIX" \
  --make-bed \
  --out "$BED_PREFIX" \
  > "${LOG_DIR}/${CHR}.plink2_make_bed.log" 2>&1

FINAL_N=$(awk 'END {print NR-1}' "${PGEN_PREFIX}.pvar")
SAMPLE_N=$(awk 'END {print NR-1}' "${PGEN_PREFIX}.psam")
{
  printf 'chrom\tpanel\tinput_variants\tafter_dr2_af\tafter_plink_qc\tsamples\n'
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$CHR" "$PANEL_UPPER" "$BEFORE_N" "$AF_DR2_N" "$FINAL_N" "$SAMPLE_N"
} > "${OUT_DIR}/qc_counts.tsv"

log "Completed ${PANEL_UPPER} QC for ${CHR}: ${FINAL_N}/${BEFORE_N} variants retained"

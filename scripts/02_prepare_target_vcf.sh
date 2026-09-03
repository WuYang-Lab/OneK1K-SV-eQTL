#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

load_config "${1:-}"
require_var HG38_FASTA
require_var CHROMOSOMES_FILE
require_command plink2
require_command bcftools
require_command tabix
require_file "$HG38_FASTA"
require_file "${HG38_FASTA}.fai"
require_file "$CHROMOSOMES_FILE"

IN_PREFIX="${PROJECT_ROOT}/results/genotypes/target.hg38"
OUT_DIR="${PROJECT_ROOT}/results/target"
WORK_DIR="${PROJECT_ROOT}/work/target_vcf"
LOG_DIR="${PROJECT_ROOT}/logs"
OUT_VCF="${OUT_DIR}/target.hg38.chr.id.vcf.gz"
mkdir -p "$OUT_DIR" "$WORK_DIR" "$LOG_DIR"
require_file "${IN_PREFIX}.bed"
require_file "${IN_PREFIX}.bim"
require_file "${IN_PREFIX}.fam"

log "Assigning GRCh38 reference alleles with PLINK 2"
plink2 \
  --bfile "$IN_PREFIX" \
  --ref-from-fa \
  --fa "$HG38_FASTA" \
  --make-pgen \
  --out "${WORK_DIR}/target.hg38.reref" \
  > "${LOG_DIR}/02_plink2_ref_from_fa.log" 2>&1

plink2 \
  --pfile "${WORK_DIR}/target.hg38.reref" \
  --export vcf bgz id-paste=iid \
  --out "${WORK_DIR}/target.hg38.reref" \
  > "${LOG_DIR}/02_plink2_export_vcf.log" 2>&1

clean_list "$CHROMOSOMES_FILE" | awk 'BEGIN{OFS="\t"} {print $1, "chr" $1}' \
  > "${WORK_DIR}/chr_map.tsv"

log "Normalizing, checking REF alleles, and indexing the target VCF"
bcftools view \
  --min-alleles 2 \
  --max-alleles 2 \
  --types snps \
  -Ou \
  "${WORK_DIR}/target.hg38.reref.vcf.gz" \
  | bcftools annotate \
      --rename-chrs "${WORK_DIR}/chr_map.tsv" \
      --set-id '%CHROM:%POS:%REF:%ALT' \
      -Ou \
  | bcftools norm \
      --fasta-ref "$HG38_FASTA" \
      --check-ref e \
      --rm-dup exact \
      -Ou \
  | bcftools sort \
      -Oz \
      -o "$OUT_VCF"

tabix -f -p vcf "$OUT_VCF"

SAMPLE_COUNT=$(bcftools query -l "$OUT_VCF" | wc -l)
VARIANT_COUNT=$(vcf_count "$OUT_VCF")
{
  printf 'metric\tvalue\n'
  printf 'samples\t%s\n' "$SAMPLE_COUNT"
  printf 'biallelic_snps\t%s\n' "$VARIANT_COUNT"
} > "${OUT_DIR}/target_vcf_summary.tsv"

log "Prepared target VCF with ${SAMPLE_COUNT} samples and ${VARIANT_COUNT} variants"
log "Output: ${OUT_VCF}"

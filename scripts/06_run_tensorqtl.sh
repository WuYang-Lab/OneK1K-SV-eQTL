#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

load_config "${1:-}"
CELL_TYPE="${2:-}"
GENOTYPE_LABEL=$(printf '%s' "${3:-sv}" | tr '[:upper:]' '[:lower:]')
MODES="${4:-cis_nominal,cis,cis_independent}"
[[ -n "$CELL_TYPE" ]] || die "Cell type is required"
require_var PHENOTYPE_TEMPLATE
require_var COVARIATE_TEMPLATE
require_command python

case "$GENOTYPE_LABEL" in
  sv)
    GENOTYPE_PREFIX="${SV_GENOTYPE_PREFIX:-${PROJECT_ROOT}/results/genotypes/sv.merged.bedfmt}"
    ;;
  vntr)
    GENOTYPE_PREFIX="${VNTR_GENOTYPE_PREFIX:-${PROJECT_ROOT}/results/genotypes/vntr.merged.bedfmt}"
    ;;
  joint)
    GENOTYPE_PREFIX="${JOINT_GENOTYPE_PREFIX:-${PROJECT_ROOT}/results/genotypes/joint/snp_sv.joint}"
    ;;
  *)
    die "Genotype label must be sv, vntr, or joint"
    ;;
esac

PHENOTYPE=$(render_cell_type_template "$PHENOTYPE_TEMPLATE" "$CELL_TYPE")
COVARIATES=$(render_cell_type_template "$COVARIATE_TEMPLATE" "$CELL_TYPE")
OUT_DIR="${PROJECT_ROOT}/results/tensorqtl/${GENOTYPE_LABEL}/${CELL_TYPE}"
PREFIX="${CELL_TYPE}.${GENOTYPE_LABEL}"
CIS_FILE="${OUT_DIR}/${PREFIX}.cis_qtl.txt.gz"
CIS_QVAL_FILE="${OUT_DIR}/${PREFIX}.cis_qtl.with_qval.txt.gz"
mkdir -p "$OUT_DIR"
require_file "$PHENOTYPE"
require_file "$COVARIATES"

python "${SCRIPT_DIR}/validate_tensorqtl_inputs.py" \
  --genotype-prefix "$GENOTYPE_PREFIX" \
  --phenotype-bed "$PHENOTYPE" \
  --covariates "$COVARIATES" \
  --output-json "${OUT_DIR}/${PREFIX}.input_validation.json"

CIS_WINDOW="${CIS_WINDOW:-1000000}"
TENSORQTL_PERMUTATIONS="${TENSORQTL_PERMUTATIONS:-10000}"
TENSORQTL_FDR="${TENSORQTL_FDR:-0.05}"
TENSORQTL_MAF="${TENSORQTL_MAF:-0.01}"
TRANS_PVAL_THRESHOLD="${TRANS_PVAL_THRESHOLD:-1e-5}"
SKIP_EXISTING="${SKIP_EXISTING:-true}"

has_mode() {
  [[ ",${MODES}," == *",$1,"* ]]
}

for requested in ${MODES//,/ }; do
  case "$requested" in
    cis_nominal|cis|cis_independent|trans) ;;
    *) die "Unsupported TensorQTL mode: ${requested}" ;;
  esac
done

if has_mode cis_nominal; then
  if [[ "$SKIP_EXISTING" == "true" ]] && compgen -G "${OUT_DIR}/${PREFIX}.cis_qtl_pairs.*.parquet" >/dev/null; then
    log "Skipping existing cis_nominal output for ${CELL_TYPE} (${GENOTYPE_LABEL})"
  else
    log "Running TensorQTL cis_nominal for ${CELL_TYPE} (${GENOTYPE_LABEL})"
    python -m tensorqtl \
      "$GENOTYPE_PREFIX" "$PHENOTYPE" "$PREFIX" \
      --covariates "$COVARIATES" \
      --mode cis_nominal \
      --window "$CIS_WINDOW" \
      --maf_threshold "$TENSORQTL_MAF" \
      --output_dir "$OUT_DIR"
  fi
fi

if has_mode cis || has_mode cis_independent; then
  if [[ -s "$CIS_FILE" && "$SKIP_EXISTING" == "true" ]]; then
    log "Skipping existing cis output for ${CELL_TYPE} (${GENOTYPE_LABEL})"
  else
    log "Running TensorQTL permutation-based cis mapping for ${CELL_TYPE} (${GENOTYPE_LABEL})"
    python -m tensorqtl \
      "$GENOTYPE_PREFIX" "$PHENOTYPE" "$PREFIX" \
      --covariates "$COVARIATES" \
      --mode cis \
      --window "$CIS_WINDOW" \
      --permutations "$TENSORQTL_PERMUTATIONS" \
      --fdr "$TENSORQTL_FDR" \
      --maf_threshold "$TENSORQTL_MAF" \
      --output_dir "$OUT_DIR"
  fi
  require_file "$CIS_FILE"

  CIS_HEADER=$(gzip -cd "$CIS_FILE" | awk 'NR == 1 {print}')
  if awk -v header="$CIS_HEADER" 'BEGIN {
    n=split(header, x, "\t"); for (i=1; i<=n; i++) if (x[i] == "qval") exit 0; exit 1
  }'; then
    CIS_FOR_INDEPENDENT="$CIS_FILE"
  else
    require_command Rscript
    log "TensorQTL cis output lacks qval; adding it with the compatibility helper"
    Rscript "${SCRIPT_DIR}/add_qvalues.R" "$CIS_FILE" "$CIS_QVAL_FILE"
    CIS_FOR_INDEPENDENT="$CIS_QVAL_FILE"
  fi
fi

if has_mode cis_independent; then
  INDEPENDENT_FILE="${OUT_DIR}/${PREFIX}.cis_independent_qtl.txt.gz"
  if [[ -s "$INDEPENDENT_FILE" && "$SKIP_EXISTING" == "true" ]]; then
    log "Skipping existing cis_independent output for ${CELL_TYPE} (${GENOTYPE_LABEL})"
  else
    log "Running TensorQTL conditionally independent cis mapping"
    python -m tensorqtl \
      "$GENOTYPE_PREFIX" "$PHENOTYPE" "$PREFIX" \
      --covariates "$COVARIATES" \
      --mode cis_independent \
      --cis_output "$CIS_FOR_INDEPENDENT" \
      --window "$CIS_WINDOW" \
      --permutations "$TENSORQTL_PERMUTATIONS" \
      --fdr "$TENSORQTL_FDR" \
      --maf_threshold "$TENSORQTL_MAF" \
      --output_dir "$OUT_DIR"
  fi
fi

if has_mode trans; then
  TRANS_FILE="${OUT_DIR}/${PREFIX}.trans_qtl_pairs.parquet"
  if [[ -s "$TRANS_FILE" && "$SKIP_EXISTING" == "true" ]]; then
    log "Skipping existing trans output for ${CELL_TYPE} (${GENOTYPE_LABEL})"
  else
    log "Running TensorQTL trans mapping"
    python -m tensorqtl \
      "$GENOTYPE_PREFIX" "$PHENOTYPE" "$PREFIX" \
      --covariates "$COVARIATES" \
      --mode trans \
      --pval_threshold "$TRANS_PVAL_THRESHOLD" \
      --maf_threshold "$TENSORQTL_MAF" \
      --output_dir "$OUT_DIR"
  fi
fi

log "TensorQTL workflow completed for ${CELL_TYPE} (${GENOTYPE_LABEL})"

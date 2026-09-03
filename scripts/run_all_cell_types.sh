#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

load_config "${1:-}"
GENOTYPE_LABEL="${2:-sv}"
MODES="${3:-cis_nominal,cis,cis_independent}"
require_var CELL_TYPES_FILE
require_file "$CELL_TYPES_FILE"

while read -r CELL_TYPE; do
  log "Starting ${CELL_TYPE}"
  bash "${SCRIPT_DIR}/06_run_tensorqtl.sh" \
    "$1" "$CELL_TYPE" "$GENOTYPE_LABEL" "$MODES"
done < <(clean_list "$CELL_TYPES_FILE")

log "Completed TensorQTL for all configured cell types"

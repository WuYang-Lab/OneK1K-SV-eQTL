#!/usr/bin/env bash

die() {
  printf '[ERROR] %s\n' "$*" >&2
  exit 1
}

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

require_file() {
  [[ -s "$1" ]] || die "Required file is missing or empty: $1"
}

require_var() {
  local name="$1"
  [[ -n "${!name:-}" ]] || die "Required configuration variable is not set: ${name}"
}

load_config() {
  local config_path="${1:-}"
  [[ -n "$config_path" ]] || die "A configuration file is required"
  require_file "$config_path"
  # shellcheck disable=SC1090
  source "$config_path"
  require_var PROJECT_ROOT
  mkdir -p "$PROJECT_ROOT"
}

clean_list() {
  awk 'NF && $1 !~ /^#/ {print $1}' "$1"
}

normalize_chr_number() {
  local value="${1#chr}"
  [[ "$value" =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || die "Expected an autosome number 1-22, received: $1"
  printf '%s\n' "$value"
}

render_cell_type_template() {
  local template="$1"
  local cell_type="$2"
  printf '%s\n' "${template//\{cell_type\}/$cell_type}"
}

vcf_count() {
  bcftools index -n "$1"
}

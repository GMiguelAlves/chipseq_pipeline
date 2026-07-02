#!/usr/bin/env bash

set -euo pipefail

CONFIG_FILE="${1:-config/pipeline_config.sh}"
if [[ "${CONFIG_FILE}" == "--config" ]]; then
  CONFIG_FILE="${2:?Usage: preflight.sh [--config FILE|FILE]}"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
load_chipseq_config "${CONFIG_FILE}"

status=0

section() {
  printf '\n== %s ==\n' "$1"
}

run_check() {
  local label="$1"
  shift
  printf '\n[%s]\n' "${label}"
  if "$@"; then
    printf 'OK: %s\n' "${label}"
  else
    printf 'FAILED: %s\n' "${label}" >&2
    status=1
  fi
}

section "Project"
printf 'Pipeline: %s\n' "${PIPELINE_NAME}"
printf 'Organism: %s\n' "${ORGANISM_NAME}"
printf 'Executor: %s\n' "${PIPELINE_EXECUTOR:-${RUN_MODE:-slurm}}"
printf 'FASTQ_DIR: %s\n' "${FASTQ_DIR}"
printf 'GENOME_FASTA: %s\n' "${GENOME_FASTA}"
printf 'ANNOTATION_FILE: %s\n' "${ANNOTATION_FILE}"
printf 'METADATA_FILE: %s\n' "${METADATA_FILE}"
printf 'WORK_ROOT: %s\n' "${WORK_ROOT}"

section "Config And Metadata"
SKIP_SLURM_CHECK=true run_check "config values and required files" bash "${SCRIPT_DIR}/validate_config.sh" "${CONFIG_FILE}"

metadata_args=(
  "${PYTHON_BIN}" "${SCRIPT_DIR}/validate_metadata.py"
  --metadata "${METADATA_FILE}"
  --fastq-dir "${FASTQ_DIR}"
  --min-replicates "${MIN_REPLICATES_DIFF}"
)
if bool_true "${REQUIRE_DIFF_REPLICATES:-false}"; then
  metadata_args+=(--require-diff-replicates)
fi
if bool_true "${ALLOW_MISSING_CONTROLS:-false}"; then
  metadata_args+=(--allow-missing-controls)
fi
run_check "metadata structure, FASTQs, controls, and replicates" "${metadata_args[@]}"

section "Experiment Summary"
sample_count="$(metadata_samples | wc -l | tr -d ' ')"
ip_count="$(metadata_ip_samples | wc -l | tr -d ' ')"
control_count="$(metadata_control_samples | wc -l | tr -d ' ')"
printf 'Samples: %s\n' "${sample_count}"
printf 'IP samples: %s\n' "${ip_count}"
printf 'Control/input samples: %s\n' "${control_count}"
printf 'Marks/factors: %s\n' "$(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="mark_or_factor") m=i; next} NR>1 && $m!=""{print $m}' "${METADATA_FILE}" | sort -u | tr '\n' ' ')"
printf 'Conditions: %s\n' "$(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="condition") c=i; next} NR>1 && $c!=""{print $c}' "${METADATA_FILE}" | sort -u | tr '\n' ' ')"

if [[ "${control_count}" -eq 0 ]]; then
  printf 'WARNING: no control/input samples are present. Use ALLOW_MISSING_CONTROLS=true only when controls are unavailable.\n' >&2
fi

section "Genome Size"
if [[ -s "${REF_DIR}/chrom.sizes" ]]; then
  printf 'MACS genome size: %s\n' "$(effective_genome_size_value "${MACS_GENOME_SIZE:-auto}" "${REF_DIR}/chrom.sizes")"
  printf 'RPGC effective genome size: %s\n' "$(effective_genome_size_value "${EFFECTIVE_GENOME_SIZE:-auto}" "${REF_DIR}/chrom.sizes")"
elif [[ -n "${EFFECTIVE_GENOME_SIZES:-}" && -s "${EFFECTIVE_GENOME_SIZES}" ]]; then
  printf 'EFFECTIVE_GENOME_SIZES is set: %s\n' "${EFFECTIVE_GENOME_SIZES}"
  printf 'The reference step will still create chrom.sizes before MACS/tracks run.\n'
else
  printf 'chrom.sizes does not exist yet; it will be created by the reference step.\n'
  printf 'For non-model organisms, review MACS_GENOME_SIZE or EFFECTIVE_GENOME_SIZES before final analysis.\n'
fi

section "Software"
run_check "required command-line tools for current config" bash "${SCRIPT_DIR}/check_install.sh" "${CONFIG_FILE}"

section "Result"
if [[ "${status}" -eq 0 ]]; then
  printf 'Preflight passed. The project is ready for dry-run or execution.\n'
else
  printf 'Preflight found issues. Fix the FAILED sections above before running the full pipeline.\n' >&2
fi

exit "${status}"

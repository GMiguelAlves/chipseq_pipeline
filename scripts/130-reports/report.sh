#!/usr/bin/env bash

set -euo pipefail

CONFIG_FILE="${1:?Usage: report.sh CONFIG}"
PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${CONFIG_FILE}")/.." && pwd)}"
SCRIPT_DIR="${PROJECT_DIR}/scripts"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
load_config "${CONFIG_FILE}"
activate_runtime

STEP="130-reports"
STEP_DIR="${REPORT_DIR:-${OUTPUT_DIR}/${STEP}}"
ensure_dir "${STEP_DIR}"

if is_done "${STEP}" "final_report"; then
  log "Final report already completed; skipping"
  exit 0
fi

if command -v "${RSCRIPT_BIN}" >/dev/null 2>&1; then
  "${RSCRIPT_BIN}" "${SCRIPT_DIR}/r/render_report.R" \
    --metadata "${METADATA_FILE}" \
    --output-dir "${WORK_ROOT:-${OUTPUT_DIR}}" \
    --report "${STEP_DIR}/chipseq_report.md"
else
  {
    echo "# ChIP-seq report"
    echo
    echo "Rscript was not available, so only a minimal report was generated."
    echo
    echo "Generated files are under: ${WORK_ROOT:-${OUTPUT_DIR}}"
  } > "${STEP_DIR}/chipseq_report.md"
fi

mark_done "${STEP}" "final_report"
log "Final report completed"

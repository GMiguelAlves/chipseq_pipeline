#!/usr/bin/env bash
set -euo pipefail

# Step-by-step runbook for PRJNA1090249 SmHP1 ChIP-seq.
# It prepares a clean run directory, downloads only the S. mansoni runs,
# writes run-specific configs, and submits the ChIP and integrative pipelines.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}"

RUN_ROOT="${RUN_ROOT:-/scratch/${USER}/smansoni_prjna1090249_smhp1}"
FASTQ_DIR="${FASTQ_DIR:-${RUN_ROOT}/data/fastq}"
THREADS="${THREADS:-8}"

GENOME_FASTA="${GENOME_FASTA:-/path/to/SM_V10/genome.fa}"
ANNOTATION_FILE="${ANNOTATION_FILE:-/path/to/SM_V10/annotation.gtf}"

METADATA_FILE="${METADATA_FILE:-${PROJECT_DIR}/config/metadata_prjna1090249_smhp1.tsv}"
CHIP_CONFIG="${CHIP_CONFIG:-${RUN_ROOT}/config/prjna1090249_smhp1.pipeline_config.sh}"

INTEGRATIVE_PROJECT_DIR="${INTEGRATIVE_PROJECT_DIR:-/path/to/integrateseq_pipeline}"
INTEGRATIVE_RUN_ROOT="${INTEGRATIVE_RUN_ROOT:-/scratch/${USER}/integrateseq_prjna1090249_smhp1}"
INTEGRATIVE_CONFIG="${INTEGRATIVE_CONFIG:-${INTEGRATIVE_RUN_ROOT}/config/prjna1090249_smhp1.integrative_config.sh}"
RNASEQ_RESULTS_DIR="${RNASEQ_RESULTS_DIR:-/path/to/rnaseq_results}"

SMANSONI_RUNS=(
  SRR28402826  # cercariae SmHP1 rep1
  SRR28402825  # cercariae SmHP1 rep2
  SRR28402822  # cercariae SmHP1 rep3
  SRR28402821  # cercariae input
  SRR28402820  # sporocysts SmHP1 rep1
  SRR28402819  # sporocysts SmHP1 rep2
  SRR28402818  # sporocysts SmHP1 rep3
  SRR28402817  # sporocysts input
)

DROSOPHILA_QC_RUNS=(
  SRR28402816
  SRR28402815
  SRR28402824
  SRR28402823
)

usage() {
  cat <<'USAGE'
Usage:
  bash examples/prjna1090249_smhp1_step_by_step.sh <step>

Required exports before running on the server:
  export GENOME_FASTA=/path/to/SM_V10/genome.fa
  export ANNOTATION_FILE=/path/to/SM_V10/annotation.gtf
  export RUN_ROOT=/scratch/$USER/smansoni_prjna1090249_smhp1

Required for the integrative step:
  export INTEGRATIVE_PROJECT_DIR=/path/to/integrateseq_pipeline
  export INTEGRATIVE_RUN_ROOT=/scratch/$USER/integrateseq_prjna1090249_smhp1
  export RNASEQ_RESULTS_DIR=/path/to/rnaseq_results

Steps:
  00-setup                  Create folders and write the ChIP config
  01-download               Download the 8 S. mansoni FASTQs from PRJNA1090249
  02-validate               Check metadata, references, FASTQs, and controls
  03-dry-run-chip           Show the full ChIP pipeline plan
  04-run-chip-all           Submit/resume the full ChIP pipeline with SLURM
  04-print-chip-steps       Print one-step-at-a-time ChIP commands
  05-check-chip             Inspect ChIP outputs and check for unknown_ChIP
  06-write-integrative      Write the integrative config pointing to this ChIP run
  07-dry-run-integrative    Show the full integrative pipeline plan
  08-run-integrative-all    Submit/resume the full integrative pipeline with SLURM
  08-print-integrative-steps Print one-step-at-a-time integrative commands
  09-check-integrative      Inspect report/figures and check for unknown_ChIP
  all-preflight             Run setup, download, validate, and dry-run checks only

Notes:
  - PRJNA1090249 also has Drosophila control/QC runs. They are intentionally
    excluded here because this S. mansoni pipeline aligns to the S. mansoni genome.
  - SmHP1 is written explicitly as mark_or_factor, so the pipeline should not
    collapse these samples to unknown_ChIP.
USAGE
}

mkdirs() {
  mkdir -p "${RUN_ROOT}/config" "${FASTQ_DIR}" "${RUN_ROOT}/sra"
}

write_chip_config() {
  mkdirs
  cat > "${RUN_ROOT}/config/disabled_user_settings.sh" <<'EOF'
# Intentionally empty. This prevents user-level defaults from changing this run.
EOF
  cat > "${RUN_ROOT}/config/disabled_pipeline_config.local.sh" <<'EOF'
# Intentionally empty. This prevents local defaults from changing this run.
EOF
  cat > "${CHIP_CONFIG}" <<EOF
#!/usr/bin/env bash
export USER_SETTINGS_FILE="${RUN_ROOT}/config/disabled_user_settings.sh"
export LOCAL_CONFIG="${RUN_ROOT}/config/disabled_pipeline_config.local.sh"

export PROJECT_DIR="${PROJECT_DIR}"
export FASTQ_DIR="${FASTQ_DIR}"
export OUTPUT_DIR="${RUN_ROOT}"
export WORK_ROOT="${RUN_ROOT}"
export METADATA_FILE="${METADATA_FILE}"
export GENOME_FASTA="${GENOME_FASTA}"
export ANNOTATION_FILE="${ANNOTATION_FILE}"

export PIPELINE_EXECUTOR="\${PIPELINE_EXECUTOR:-slurm}"
export RUN_MODE="\${RUN_MODE:-\${PIPELINE_EXECUTOR}}"
export PEAK_CALLER="\${PEAK_CALLER:-macs3}"
export PEAK_TYPE="\${PEAK_TYPE:-narrow}"
export ALLOW_MISSING_CONTROLS=false
export GROUP_COLUMNS="condition,mark_or_factor"
export MIN_REPLICATES_DIFF=2
export REQUIRE_DIFF_REPLICATES=false

source "${PROJECT_DIR}/config/pipeline_config.sh"
EOF
  chmod +x "${CHIP_CONFIG}"
  echo "Wrote ${CHIP_CONFIG}"
}

download_fastqs() {
  mkdirs
  command -v prefetch >/dev/null 2>&1 || { echo "Missing prefetch from SRA Toolkit"; exit 1; }
  command -v fasterq-dump >/dev/null 2>&1 || { echo "Missing fasterq-dump from SRA Toolkit"; exit 1; }

  for run in "${SMANSONI_RUNS[@]}"; do
    if [[ -s "${FASTQ_DIR}/${run}.fastq.gz" ]]; then
      echo "Already present: ${FASTQ_DIR}/${run}.fastq.gz"
      continue
    fi

    echo "Downloading ${run}"
    prefetch "${run}" --output-directory "${RUN_ROOT}/sra"
    fasterq-dump "${run}" --outdir "${FASTQ_DIR}" --threads "${THREADS}" --split-files --skip-technical

    if [[ -s "${FASTQ_DIR}/${run}_1.fastq" && ! -s "${FASTQ_DIR}/${run}_2.fastq" ]]; then
      mv "${FASTQ_DIR}/${run}_1.fastq" "${FASTQ_DIR}/${run}.fastq"
    fi

    if [[ ! -s "${FASTQ_DIR}/${run}.fastq" ]]; then
      echo "Expected single-end FASTQ not found for ${run}" >&2
      exit 1
    fi

    if command -v pigz >/dev/null 2>&1; then
      pigz -p "${THREADS}" "${FASTQ_DIR}/${run}.fastq"
    else
      gzip "${FASTQ_DIR}/${run}.fastq"
    fi
  done
}

validate_inputs() {
  [[ -s "${METADATA_FILE}" ]] || { echo "Missing metadata: ${METADATA_FILE}" >&2; exit 1; }
  [[ -s "${GENOME_FASTA}" ]] || { echo "Missing genome FASTA: ${GENOME_FASTA}" >&2; exit 1; }
  [[ -s "${ANNOTATION_FILE}" ]] || { echo "Missing annotation GTF: ${ANNOTATION_FILE}" >&2; exit 1; }

  awk -F'\t' 'NR > 1 {print $2}' "${METADATA_FILE}" | while read -r fastq; do
    [[ -s "${FASTQ_DIR}/${fastq}" ]] || { echo "Missing FASTQ: ${FASTQ_DIR}/${fastq}" >&2; exit 1; }
  done

  python3 "${PROJECT_DIR}/scripts/validate_metadata.py" \
    --metadata "${METADATA_FILE}" \
    --fastq-dir "${FASTQ_DIR}" \
    --min-replicates 2
}

dry_run_chip() {
  write_chip_config
  bash "${PROJECT_DIR}/chipseq_pipeline.sh" --config "${CHIP_CONFIG}" --all --slurm --dry-run
}

run_chip_all() {
  write_chip_config
  bash "${PROJECT_DIR}/chipseq_pipeline.sh" --config "${CHIP_CONFIG}" --all --slurm --resume
}

print_chip_steps() {
  write_chip_config
  for step in reference qc trim align filter bam_qc peaks consensus differential annotate tracks report; do
    echo "bash ${PROJECT_DIR}/chipseq_pipeline.sh --config ${CHIP_CONFIG} --step ${step} --slurm --resume"
  done
  echo "Run each command only after the previous SLURM jobs have completed successfully."
}

check_chip() {
  echo "Metadata factors:"
  awk -F'\t' 'NR > 1 {print $6}' "${METADATA_FILE}" | sort -u

  echo
  echo "Annotated peak files:"
  find "${RUN_ROOT}/090-peak-annotation" -name '*.annotated.tsv*' -print 2>/dev/null | sort || true

  echo
  echo "Consensus peak outputs:"
  find "${RUN_ROOT}/110-consensus-peaks" -type f -print 2>/dev/null | sort | head -100 || true

  echo
  echo "unknown_ChIP search:"
  grep -R "unknown_ChIP" \
    "${RUN_ROOT}/090-peak-annotation" \
    "${RUN_ROOT}/110-consensus-peaks" \
    "${RUN_ROOT}/120-differential-binding" 2>/dev/null || true
}

write_integrative_config() {
  mkdir -p "${INTEGRATIVE_RUN_ROOT}/config"
  cat > "${INTEGRATIVE_CONFIG}" <<EOF
#!/usr/bin/env bash
export INTEGRATION_OUTPUT_DIR="${INTEGRATIVE_RUN_ROOT}"
export RNASEQ_RESULTS_DIR="${RNASEQ_RESULTS_DIR}"
export CHIPSEQ_RESULTS_DIR="${RUN_ROOT}"
export ANNOTATION_FILE="${ANNOTATION_FILE}"

export CHIP_METADATA_FILE="${METADATA_FILE}"
export CHIP_ANNOTATED_PEAKS_GLOB="${RUN_ROOT}/090-peak-annotation/*.annotated.tsv*"
export CHIP_PEAK_BED_GLOB="${RUN_ROOT}/110-consensus-peaks/groups/*.consensus.bed*"
export CHIP_PEAK_COUNT_GLOB="${RUN_ROOT}/110-consensus-peaks/counts/*.counts.tsv*"
export CHIP_DIFF_BINDING_FILE="${RUN_ROOT}/120-differential-binding/differential_binding_results.tsv.gz"

source "${INTEGRATIVE_PROJECT_DIR}/config/pipeline_config.sh"

# Re-apply run-specific paths after sourcing defaults/local settings.
export INTEGRATION_OUTPUT_DIR="${INTEGRATIVE_RUN_ROOT}"
export RNASEQ_RESULTS_DIR="${RNASEQ_RESULTS_DIR}"
export CHIPSEQ_RESULTS_DIR="${RUN_ROOT}"
export ANNOTATION_FILE="${ANNOTATION_FILE}"
export CHIP_METADATA_FILE="${METADATA_FILE}"
export CHIP_ANNOTATED_PEAKS_GLOB="${RUN_ROOT}/090-peak-annotation/*.annotated.tsv*"
export CHIP_PEAK_BED_GLOB="${RUN_ROOT}/110-consensus-peaks/groups/*.consensus.bed*"
export CHIP_PEAK_COUNT_GLOB="${RUN_ROOT}/110-consensus-peaks/counts/*.counts.tsv*"
if [[ -s "${RUN_ROOT}/120-differential-binding/differential_binding_results.tsv.gz" ]]; then
  export CHIP_DIFF_BINDING_FILE="${RUN_ROOT}/120-differential-binding/differential_binding_results.tsv.gz"
else
  export CHIP_DIFF_BINDING_FILE="${RUN_ROOT}/120-differential-binding/differential_binding_results.tsv"
fi
EOF
  chmod +x "${INTEGRATIVE_CONFIG}"
  echo "Wrote ${INTEGRATIVE_CONFIG}"
}

dry_run_integrative() {
  write_integrative_config
  bash "${INTEGRATIVE_PROJECT_DIR}/integrative_pipeline.sh" --config "${INTEGRATIVE_CONFIG}" --all --dry-run
}

run_integrative_all() {
  write_integrative_config
  bash "${INTEGRATIVE_PROJECT_DIR}/integrative_pipeline.sh" --config "${INTEGRATIVE_CONFIG}" --all --slurm --resume
}

print_integrative_steps() {
  write_integrative_config
  for step in validate prepare harmonize map-peaks summarize-rna summarize-chip integrate score visualize functional report; do
    echo "bash ${INTEGRATIVE_PROJECT_DIR}/integrative_pipeline.sh --config ${INTEGRATIVE_CONFIG} --step ${step} --slurm --resume"
  done
  echo "Run each command only after the previous SLURM jobs have completed successfully."
}

check_integrative() {
  echo "Reports:"
  find "${INTEGRATIVE_RUN_ROOT}" -maxdepth 3 -type f \( -name '*.html' -o -name '*.pdf' \) -print 2>/dev/null | sort || true

  echo
  echo "Visualizations:"
  find "${INTEGRATIVE_RUN_ROOT}" -maxdepth 3 -type f \( -name '*.png' -o -name '*.svg' \) -print 2>/dev/null | sort | head -100 || true

  echo
  echo "unknown_ChIP search:"
  grep -R "unknown_ChIP" "${INTEGRATIVE_RUN_ROOT}" 2>/dev/null || true
}

case "${1:-help}" in
  00-setup)
    write_chip_config
    ;;
  01-download)
    download_fastqs
    ;;
  02-validate)
    validate_inputs
    ;;
  03-dry-run-chip)
    dry_run_chip
    ;;
  04-run-chip-all)
    run_chip_all
    ;;
  04-print-chip-steps)
    print_chip_steps
    ;;
  05-check-chip)
    check_chip
    ;;
  06-write-integrative)
    write_integrative_config
    ;;
  07-dry-run-integrative)
    dry_run_integrative
    ;;
  08-run-integrative-all)
    run_integrative_all
    ;;
  08-print-integrative-steps)
    print_integrative_steps
    ;;
  09-check-integrative)
    check_integrative
    ;;
  all-preflight)
    write_chip_config
    download_fastqs
    validate_inputs
    dry_run_chip
    write_integrative_config
    dry_run_integrative
    ;;
  help|-h|--help)
    usage
    ;;
  *)
    usage
    echo "Unknown step: ${1}" >&2
    exit 2
    ;;
esac

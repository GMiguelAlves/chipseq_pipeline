#!/usr/bin/env bash
# Central configuration for the generic ChIP-seq pipeline.
# Copy this file or edit it in place before running chipseq_pipeline.sh.

set -euo pipefail

# Project layout
PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# Optional user-local overrides. This file is ignored by Git so users can keep
# HPC paths, accounts, and reference locations outside version control.
LOCAL_CONFIG="${LOCAL_CONFIG:-${PROJECT_DIR}/config/pipeline_config.local.sh}"
if [[ -f "${LOCAL_CONFIG}" ]]; then
  # shellcheck source=/dev/null
  source "${LOCAL_CONFIG}"
fi

FASTQ_DIR="${FASTQ_DIR:-${PROJECT_DIR}/data/fastq}"
REFERENCE_DIR="${REFERENCE_DIR:-${PROJECT_DIR}/reference}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}}"
METADATA_FILE="${METADATA_FILE:-${PROJECT_DIR}/config/metadata.tsv}"

# Organism-dependent reference files. These must be supplied by the user.
GENOME_FASTA="${GENOME_FASTA:-${REFERENCE_DIR}/genome.fa}"
ANNOTATION_FILE="${ANNOTATION_FILE:-${REFERENCE_DIR}/annotation.gtf}"
BLACKLIST_BED="${BLACKLIST_BED:-}"
EFFECTIVE_GENOME_SIZES="${EFFECTIVE_GENOME_SIZES:-}"
FUNCTIONAL_ANNOTATION="${FUNCTIONAL_ANNOTATION:-}"

# Runtime
RUN_MODE="${RUN_MODE:-slurm}"              # slurm or local
THREADS="${THREADS:-8}"
DEFAULT_MEM="${DEFAULT_MEM:-32G}"
SLURM_TIME="${SLURM_TIME:-12:00:00}"
SLURM_PARTITION="${SLURM_PARTITION:-compute}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"

# Software environment. ENV_BACKEND can be none, conda, or apptainer.
ENV_BACKEND="${ENV_BACKEND:-conda}"
CONDA_ENV="${CONDA_ENV:-chipseq}"
CONTAINER_IMAGE="${CONTAINER_IMAGE:-}"

# FASTQ and trimming
READ_LAYOUT="${READ_LAYOUT:-metadata}"     # metadata, paired, or single
TRIM_TOOL="${TRIM_TOOL:-fastp}"            # fastp or trim_galore
FASTP_OPTS="${FASTP_OPTS:---detect_adapter_for_pe --qualified_quality_phred 20 --length_required 20}"
TRIM_GALORE_OPTS="${TRIM_GALORE_OPTS:---quality 20 --length 20}"

# Alignment
ALIGNER="${ALIGNER:-bowtie2}"              # bowtie2 or bwa
BOWTIE2_INDEX_PREFIX="${BOWTIE2_INDEX_PREFIX:-${OUTPUT_DIR}/010-reference/bowtie2/genome}"
BWA_INDEX_PREFIX="${BWA_INDEX_PREFIX:-${GENOME_FASTA}}"
BOWTIE2_BUILD_OPTS="${BOWTIE2_BUILD_OPTS:-}"
BOWTIE2_OPTS="${BOWTIE2_OPTS:---very-sensitive}"
BWA_OPTS="${BWA_OPTS:-}"

# BAM filtering
MIN_MAPQ="${MIN_MAPQ:-30}"
REMOVE_SECONDARY_SUPPLEMENTARY="${REMOVE_SECONDARY_SUPPLEMENTARY:-true}"
REMOVE_DUPLICATES="${REMOVE_DUPLICATES:-true}"
DEDUP_TOOL="${DEDUP_TOOL:-samtools}"       # samtools or picard
PICARD_CMD="${PICARD_CMD:-picard}"

# Peak calling
PEAK_CALLER="${PEAK_CALLER:-macs3}"        # macs3 or macs2
PEAK_TYPE="${PEAK_TYPE:-auto}"             # auto, narrow, or broad
BROAD_MARK_REGEX="${BROAD_MARK_REGEX:-H3K27me3|H3K36me3|H3K9me3|H3K79me2}"
MACS_QVALUE="${MACS_QVALUE:-0.01}"
MACS_EXTRA_OPTS="${MACS_EXTRA_OPTS:-}"
MACS_GENOME_SIZE="${MACS_GENOME_SIZE:-auto}"

# Annotation
PROMOTER_UPSTREAM="${PROMOTER_UPSTREAM:-2000}"
PROMOTER_DOWNSTREAM="${PROMOTER_DOWNSTREAM:-500}"

# Tracks
BIGWIG_NORMALIZATION="${BIGWIG_NORMALIZATION:-CPM}"  # CPM or RPGC
EFFECTIVE_GENOME_SIZE="${EFFECTIVE_GENOME_SIZE:-auto}"
BIN_SIZE="${BIN_SIZE:-10}"

# Differential binding
GROUP_COLUMNS="${GROUP_COLUMNS:-condition,mark_or_factor}"
DIFF_CONTRASTS="${DIFF_CONTRASTS:-}"        # Example: treated:control,drug:vehicle
MIN_REPLICATES_DIFF="${MIN_REPLICATES_DIFF:-2}"
REQUIRE_DIFF_REPLICATES="${REQUIRE_DIFF_REPLICATES:-false}"

# Output safety
OVERWRITE="${OVERWRITE:-false}"
CREATE_DONE_FILES="${CREATE_DONE_FILES:-true}"

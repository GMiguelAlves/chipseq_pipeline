#!/usr/bin/env bash

# Copy this file to config/user_settings.sh and edit only this small block.
# Leave config/pipeline_config.sh alone unless you need advanced behavior.

# 1) Name your analysis.
export PIPELINE_NAME="chipseq_project"
export ORGANISM_NAME="My organism"

# 2) Point to input data and metadata.
# FASTQ paths in config/metadata.tsv may be absolute or relative to FASTQ_DIR.
export FASTQ_DIR="/path/to/fastq"
export METADATA_FILE="${PROJECT_DIR}/config/metadata.tsv"

# 3) Provide reference files for the organism.
export REFERENCE_DIR="/path/to/reference"
export GENOME_FASTA="${REFERENCE_DIR}/genome.fa"
export ANNOTATION_FILE="${REFERENCE_DIR}/annotation.gtf"

# Optional organism/project annotations.
export BLACKLIST_BED=""
export EFFECTIVE_GENOME_SIZES=""
export FUNCTIONAL_ANNOTATION=""

# 4) Choose where outputs are written.
# Keeping OUTPUT_DIR as PROJECT_DIR reproduces the numbered directory layout.
export OUTPUT_DIR="${PROJECT_DIR}"

# Large TSV-like results are written as .tsv.gz by default, like the RNA-seq
# pipeline. Set this to 0 only if you need plain text tables.
export PIPELINE_COMPRESS_RESULTS=1

# 5) Point to Conda on the Slurm server.
# Absolute paths are safest. Relative paths are resolved from this file's
# directory, so "../miniconda3" means "<project>/miniconda3".
export CONDA_BASE="/path/to/miniconda3"

# 6) Choose tools.
export ALIGNER="bowtie2"          # bowtie2 or bwa
export TRIM_TOOL="fastp"          # fastp or trim_galore
export PEAK_CALLER="macs3"        # macs3 or macs2
export PEAK_TYPE="auto"           # auto, narrow, or broad
# Set to true only when matched input/control FASTQs are not available.
# MACS will then run without -c for IP samples with an empty control_id.
export ALLOW_MISSING_CONTROLS="false"

# 7) Choose where jobs run.
export PIPELINE_EXECUTOR="slurm"  # Use "local" to run without sbatch.
export LOCAL_CPUS_PER_TASK=8      # Used only when PIPELINE_EXECUTOR="local".

# 8) Slurm resources.
export THREADS=8
export MEMORY="32G"
export SLURM_TIME="12:00:00"
export SLURM_PARTITION="compute"
export SLURM_ACCOUNT=""
export SLURM_QOS=""

# 9) Conda environment names. Change only if your server uses other names.
export ENV_BACKEND="conda"
export CHIPSEQ_TOOLS_ENV="chipseq"
export PYTHON_ENV="chipseq"
export R_ANALYSIS_ENV="chipseq"
export CONDA_ENV="${CHIPSEQ_TOOLS_ENV}"

# 10) ChIP-seq parameters.
export READ_LAYOUT="metadata"     # metadata, paired, or single
export MIN_MAPQ=30
export REMOVE_DUPLICATES="true"
export PROMOTER_UPSTREAM=2000
export PROMOTER_DOWNSTREAM=500
export BIGWIG_NORMALIZATION="CPM"
export DIFF_CONTRASTS=""          # Example: treated:control,drug:vehicle

# 11) Storage mode.
# full: keep everything. Best for debugging and reruns.
# balanced: remove temporary uncompressed reference copies and individual
#           FastQC folders after the final report.
# minimal: also remove trimmed FASTQs after final BAMs/reports exist.
export PIPELINE_STORAGE_MODE="full"

# Generic ChIP-seq Pipeline for HPC/Slurm

This repository contains a modular ChIP-seq pipeline designed for reproducible
analysis on HPC systems. It is organism agnostic: genome FASTA, GTF/GFF3,
blacklist, genome sizes, metadata, and all organism-specific choices are
configured by the user.

The pipeline starts from raw FASTQ files and produces QC reports, trimmed reads,
aligned and filtered BAMs, MACS peaks, consensus peak sets, gzipped count matrices,
differential binding outputs, peak annotations, bigWig tracks, and a final
Markdown report.

## Quick Start After Cloning

Clone the repository on the HPC login node or in a shared project directory:

```bash
git clone <repository-url> chipseq-pipeline
cd chipseq-pipeline
```

Create the suggested Conda environment:

```bash
conda env create -f envs/chipseq.yml
conda activate chipseq
```

Create local, untracked configuration files:

```bash
bash scripts/init_project.sh
```

Edit `config/user_settings.sh` with HPC paths, reference files, Slurm
settings, and software environment. Edit `config/metadata.tsv` with the samples
for the experiment.

Check the installed tools for the current config:

```bash
bash scripts/check_install.sh
```

Validate the pipeline graph before submitting jobs:

```bash
bash chipseq_pipeline.sh --all --dry-run
```

Submit the complete workflow:

```bash
bash chipseq_pipeline.sh --all
```

## Directory Layout

- `config/`: central configuration and metadata template
- `scripts/`: independent step scripts and shared functions
- `scripts/010-reference`, `scripts/030-qc-fastq`, ...: RNA-seq-style step wrappers
- `scripts/r/`: R scripts for annotation, differential analysis, and reporting
- `slurm/`: Slurm notes
- `envs/`: Conda environment template
- `000-logs/`: Slurm stdout/stderr logs by step and sample
- `010-reference/`: genome indexes, chromosome sizes, annotation BEDs
- `030-qc-fastq/`: FastQC and MultiQC outputs
- `040-trimming/`: trimmed FASTQs and trimming reports
- `050-alignment/`: sorted BAMs from Bowtie2 or BWA
- `060-filtering/`: final filtered BAMs and BAM metrics
- `070-qc-alignment/`: alignment QC outputs
- `080-peak-calling/`: MACS2/MACS3 peak calls
- `090-peak-annotation/`: annotated peak tables
- `100-tracks/`: bigWig tracks
- `110-consensus-peaks/`: consensus BEDs and count matrices
- `120-differential-binding/`: differential enrichment tables and plots
- `130-reports/`: final report

## Requirements

Install the suggested Conda environment:

```bash
conda env create -f envs/chipseq.yml
conda activate chipseq
```

Core tools: gzip, FastQC, MultiQC, fastp or Trim Galore, Bowtie2 or BWA, samtools,
bedtools, deepTools, MACS2/MACS3, Python 3, and R. DESeq2 is used when
available for differential binding; otherwise the R script writes an
exploratory fallback table.

## Configure

The versioned default config is `config/pipeline_config.sh`. It works like the
RNA-seq pipeline: most users should leave it alone and put local paths and HPC
settings in:

```text
config/user_settings.sh
```

That file is ignored by Git and is sourced automatically by
`config/pipeline_config.sh`.

Required organism-specific inputs:

- `FASTQ_DIR`
- `GENOME_FASTA`
- `ANNOTATION_FILE`
- `METADATA_FILE`

`GENOME_FASTA` and `ANNOTATION_FILE` may be plain files or gzip files such as
`.fa.gz`, `.gtf.gz`, or `.gff3.gz`. Gzipped references are decompressed into
`010-reference/data/` before indexing.

Optional inputs:

- `BLACKLIST_BED`
- `EFFECTIVE_GENOME_SIZES`
- `FUNCTIONAL_ANNOTATION`

Important execution settings:

- `PIPELINE_EXECUTOR=slurm` or `PIPELINE_EXECUTOR=local`
- `THREADS`, `MEMORY`, `SLURM_TIME`, `SLURM_PARTITION`, `SLURM_ACCOUNT`
- `ENV_BACKEND=conda`, `none`, or `apptainer`
- `ALIGNER=bowtie2` or `bwa`
- `TRIM_TOOL=fastp` or `trim_galore`
- `PEAK_CALLER=macs3` or `macs2`
- `PEAK_TYPE=auto`, `narrow`, or `broad`
- `PROMOTER_UPSTREAM` and `PROMOTER_DOWNSTREAM`
- `PIPELINE_COMPRESS_RESULTS=1` to write large TSV-like outputs as `.tsv.gz`
- `PIPELINE_STORAGE_MODE=full|balanced|minimal`

## Metadata

Copy the versioned template and edit the untracked metadata file:

```bash
cp config/metadata_template.tsv config/metadata.tsv
```

Required columns:

```text
sample_id fastq_1 fastq_2 layout assay mark_or_factor condition replicate batch treatment control_id is_control organism genome_id
```

Rules enforced before execution:

- `sample_id` must be unique.
- `layout` must be `single` or `paired`.
- If `layout=paired`, `fastq_2` is required.
- If `layout=single`, `fastq_2` may be empty.
- If `is_control=false`, `control_id` must point to an existing control sample.
- If `is_control=true`, the sample is treated as input/control and peak calling is skipped.
- Replicate counts are checked for differential binding readiness.

## Run

Complete run with Slurm dependencies:

```bash
bash chipseq_pipeline.sh --all
```

Inspect the job graph without submitting:

```bash
bash chipseq_pipeline.sh --all --dry-run
```

Run locally:

```bash
bash chipseq_pipeline.sh --all --mode local
```

Preferred RNA-seq-style local shortcut:

```bash
bash chipseq_pipeline.sh --all --local
```

Run one or more steps:

```bash
bash chipseq_pipeline.sh --step reference
bash chipseq_pipeline.sh --step qc --step trim
bash chipseq_pipeline.sh --step peaks --step consensus --step annotate
```

Supported steps:

```text
reference qc trim align filter bam_qc peaks consensus differential annotate tracks report cleanup
```

By default the pipeline resumes from `.done` files and does not overwrite
completed steps. Use `--force` to rerun selected steps.

## Workflow

Full execution order:

```text
reference
qc -> trim -> post-trim qc
trim + reference -> align -> filter -> bam_qc
filter -> peaks -> consensus -> differential
peaks + consensus + reference -> annotate
filter -> tracks
all aggregate outputs -> report
```

In Slurm mode, sample-level jobs are submitted independently and downstream
steps use `--dependency=afterok`. The orchestrator submits jobs in the same
style as the RNA-seq pipeline, using `sbatch --chdir`, `--export`,
`--parsable`, and a shared `submit_sbatch` helper.

## Outputs

Key files:

- `010-reference/chrom.sizes`
- `010-reference/annotation/{genes,exons,introns,promoters,downstream}.bed`
- `050-alignment/<sample>/<sample>.sorted.bam`
- `060-filtering/<sample>/<sample>.filtered.bam`
- `080-peak-calling/<sample>/*_peaks.narrowPeak` or `*_peaks.broadPeak`
- `110-consensus-peaks/groups/*.consensus.bed`
- `110-consensus-peaks/counts/*.counts.tsv.gz`
- `090-peak-annotation/*.annotated.tsv.gz`
- `100-tracks/*.bw`
- `120-differential-binding/differential_binding_results.tsv.gz`
- `130-reports/chipseq_report.md`

When `PIPELINE_COMPRESS_RESULTS=0`, TSV-like outputs use plain `.tsv`.

## Storage Modes

- `full`: keep every generated file.
- `balanced`: after the final report, remove temporary uncompressed reference
  copies made from `.gz` inputs and individual FastQC folders.
- `minimal`: additionally remove trimmed FASTQs after downstream BAMs and
  reports exist.

Run cleanup manually with:

```bash
bash chipseq_pipeline.sh --step cleanup
```

## Recovering Failed Jobs

1. Inspect `000-logs/<step>/<sample>.err`.
2. Fix the cause, such as a missing FASTQ, missing tool, bad reference path, or
   insufficient memory/time.
3. Rerun the failed step:

```bash
bash chipseq_pipeline.sh --step <step>
```

Completed outputs with `.done` files are skipped. To rerun after deleting or
replacing outputs, use:

```bash
bash chipseq_pipeline.sh --step <step> --force
```

## ChIP-seq Practices

- Use matched input/control samples for each IP whenever possible.
- Keep biological replicates for each condition and mark/factor.
- Check raw and post-trim FastQC before interpreting peaks.
- Inspect mapping rate, duplication, blacklist fraction, fragment size, and FRiP.
- Prefer narrow peaks for transcription factors and punctate marks.
- Prefer broad peaks for broad histone marks such as H3K27me3.
- Review bigWig tracks in IGV or UCSC before final biological interpretation.

## Limitations

- The pipeline expects the user to provide compatible FASTA and GTF/GFF3 files.
- Differential binding requires adequate replicates and a meaningful contrast
  design in metadata.
- The base R annotation is intentionally portable; large genomes may benefit
  from a GenomicRanges or ChIPseeker-based annotation workflow.
- Apptainer/Singularity support is command-wrapper based and assumes reference
  paths are visible inside the container.

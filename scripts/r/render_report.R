#!/usr/bin/env Rscript

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    value <- args[[i + 1]]
    out[[key]] <- value
    i <- i + 2
  }
  out
}

read_tsv_safe <- function(path, header = TRUE) {
  candidates <- unique(c(path, paste0(path, ".gz"), sub("\\.gz$", "", path)))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0 || file.info(existing[[1]])$size == 0) return(data.frame())
  con <- if (grepl("\\.gz$", existing[[1]])) gzfile(existing[[1]], "rt") else file(existing[[1]], "rt")
  on.exit(close(con), add = TRUE)
  read.table(con, sep = "\t", header = header, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
}

extract_flagstat <- function(path) {
  if (!file.exists(path)) return(c(total = NA, mapped_pct = NA, duplicate_pct = NA))
  lines <- readLines(path, warn = FALSE)
  total <- suppressWarnings(as.numeric(sub(" .*", "", lines[1])))
  mapped_line <- grep(" mapped \\(", lines, value = TRUE)
  mapped_pct <- if (length(mapped_line)) suppressWarnings(as.numeric(sub(".*\\(([^%]+)%.*", "\\1", mapped_line[1]))) else NA
  dup_line <- grep(" duplicates$", lines, value = TRUE)
  duplicate_pct <- if (length(dup_line) && !is.na(total) && total > 0) {
    suppressWarnings(as.numeric(sub(" .*", "", dup_line[1]))) / total * 100
  } else {
    NA
  }
  c(total = total, mapped_pct = mapped_pct, duplicate_pct = duplicate_pct)
}

args <- parse_args()
metadata <- read_tsv_safe(args[["metadata"]])
outdir <- args[["output-dir"]]
report <- args[["report"]]
dir.create(dirname(report), recursive = TRUE, showWarnings = FALSE)

is_control <- tolower(as.character(metadata$is_control)) %in% c("true", "1", "yes", "y")
sample_ids <- metadata$sample_id

sample_summary <- data.frame(
  sample_id = sample_ids,
  condition = metadata$condition,
  mark_or_factor = metadata$mark_or_factor,
  is_control = is_control,
  filtered_bam = file.exists(file.path(outdir, "060-filtering", sample_ids, paste0(sample_ids, ".filtered.bam"))),
  stringsAsFactors = FALSE
)

metrics <- t(vapply(sample_ids, function(sid) {
  extract_flagstat(file.path(outdir, "060-filtering", sid, paste0(sid, ".flagstat.txt")))
}, numeric(3)))
sample_summary$total_filtered_reads <- metrics[, "total"]
sample_summary$mapped_pct <- metrics[, "mapped_pct"]
sample_summary$duplicate_pct <- metrics[, "duplicate_pct"]

peak_counts <- vapply(sample_ids, function(sid) {
  files <- list.files(file.path(outdir, "080-peak-calling", sid), pattern = "_peaks\\.(narrowPeak|broadPeak)$", full.names = TRUE)
  if (length(files) == 0 || !file.exists(files[[1]])) return(NA_integer_)
  length(readLines(files[[1]], warn = FALSE))
}, integer(1))
sample_summary$peak_count <- peak_counts

annotation_summary <- read_tsv_safe(file.path(outdir, "090-peak-annotation", "annotation_summary.tsv"))
diff_results <- read_tsv_safe(file.path(outdir, "120-differential-binding", "differential_binding_results.tsv"))

recommendations <- character()
low_mapping <- sample_summary$sample_id[!is.na(sample_summary$mapped_pct) & sample_summary$mapped_pct < 70]
if (length(low_mapping) > 0) {
  recommendations <- c(recommendations, paste("Inspect low mapping samples:", paste(low_mapping, collapse = ", ")))
}
high_dup <- sample_summary$sample_id[!is.na(sample_summary$duplicate_pct) & sample_summary$duplicate_pct > 50]
if (length(high_dup) > 0) {
  recommendations <- c(recommendations, paste("High duplicate fraction detected:", paste(high_dup, collapse = ", ")))
}
zero_peaks <- sample_summary$sample_id[!is.na(sample_summary$peak_count) & sample_summary$peak_count == 0]
if (length(zero_peaks) > 0) {
  recommendations <- c(recommendations, paste("No peaks detected for:", paste(zero_peaks, collapse = ", ")))
}
if (length(recommendations) == 0) {
  recommendations <- "No automatic quality warnings were triggered."
}

con <- file(report, open = "w", encoding = "UTF-8")
writeLines(c(
  "# ChIP-seq final report",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Samples",
  "",
  paste("Total samples:", nrow(metadata)),
  paste("IP samples:", sum(!is_control)),
  paste("Control samples:", sum(is_control)),
  "",
  "## Sample metrics",
  ""
), con)
utils::write.table(sample_summary, con, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c("", "## Peak annotation", ""), con)
if (nrow(annotation_summary) > 0) {
  utils::write.table(annotation_summary, con, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  writeLines("Peak annotation summary was not found.", con)
}
writeLines(c("", "## Differential binding", ""), con)
if (nrow(diff_results) > 0) {
  writeLines(paste("Differential regions reported:", nrow(diff_results)), con)
} else {
  writeLines("No differential binding results were generated.", con)
}
writeLines(c("", "## Main output directories", "",
             "- 030-qc-fastq: raw and post-trim FastQC/MultiQC",
             "- 050-alignment: sorted BAMs from the aligner",
             "- 060-filtering: final filtered BAMs and BAM metrics",
             "- 080-peak-calling: MACS peak calls",
             "- 090-peak-annotation: annotated peak tables",
             "- 100-tracks: bigWig tracks",
             "- 110-consensus-peaks: merged peaks and count matrices",
             "- 120-differential-binding: differential enrichment tables and plots",
             "",
             "## Automatic recommendations",
             ""), con)
writeLines(paste("-", recommendations), con)
close(con)

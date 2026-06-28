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

first_existing <- function(path) {
  candidates <- unique(c(path, paste0(path, ".gz"), sub("\\.gz$", "", path)))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) return(NA_character_)
  existing[[1]]
}

open_table <- function(path, mode = "rt") {
  if (grepl("\\.gz$", path)) gzfile(path, mode) else file(path, mode)
}

read_tsv_safe <- function(path, header = TRUE) {
  existing <- first_existing(path)
  if (is.na(existing) || file.info(existing)$size == 0) return(data.frame())
  con <- open_table(existing, "rt")
  on.exit(close(con), add = TRUE)
  read.table(con, sep = "\t", header = header, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv_safe <- function(x, path) {
  con <- open_table(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE)
}

read_lines_safe <- function(path) {
  existing <- first_existing(path)
  if (is.na(existing) || file.info(existing)$size == 0) return(character())
  con <- if (grepl("\\.gz$", existing)) gzfile(existing, "rt") else file(existing, "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

count_lines_safe <- function(path) {
  existing <- first_existing(path)
  if (is.na(existing) || file.info(existing)$size == 0) return(NA_integer_)
  con <- if (grepl("\\.gz$", existing)) gzfile(existing, "rt") else file(existing, "rt")
  on.exit(close(con), add = TRUE)
  length(readLines(con, warn = FALSE))
}

num_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(length(x) == 0, NA_real_, x)
}

percent <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", paste0(format(round(x, digits), nsmall = digits, trim = TRUE), "%"))
}

fmt_num <- function(x) {
  ifelse(is.na(x), "NA", format(x, big.mark = ",", scientific = FALSE, trim = TRUE))
}

extract_flagstat <- function(path) {
  lines <- read_lines_safe(path)
  if (length(lines) == 0) {
    return(c(total = NA, mapped = NA, mapped_pct = NA, duplicates = NA, duplicate_pct = NA))
  }
  total <- suppressWarnings(as.numeric(sub(" .*", "", lines[1])))
  mapped_line <- grep(" mapped \\(", lines, value = TRUE)
  mapped <- if (length(mapped_line)) suppressWarnings(as.numeric(sub(" .*", "", mapped_line[1]))) else NA_real_
  mapped_pct <- if (length(mapped_line)) suppressWarnings(as.numeric(sub(".*\\(([^%]+)%.*", "\\1", mapped_line[1]))) else NA_real_
  dup_line <- grep(" duplicates$", lines, value = TRUE)
  duplicates <- if (length(dup_line)) suppressWarnings(as.numeric(sub(" .*", "", dup_line[1]))) else NA_real_
  duplicate_pct <- if (!is.na(duplicates) && !is.na(total) && total > 0) duplicates / total * 100 else NA_real_
  c(total = total, mapped = mapped, mapped_pct = mapped_pct, duplicates = duplicates, duplicate_pct = duplicate_pct)
}

extract_samtools_stats <- function(path) {
  lines <- read_lines_safe(path)
  keys <- c(
    total_sequences = "raw total sequences",
    filtered_sequences = "filtered sequences",
    reads_mapped = "reads mapped",
    reads_duplicated = "reads duplicated",
    reads_mq0 = "reads MQ0",
    average_length = "average length",
    insert_size_average = "insert size average"
  )
  out <- setNames(rep(NA_real_, length(keys)), names(keys))
  sn <- lines[grepl("^SN\t", lines)]
  for (nm in names(keys)) {
    hit <- grep(paste0("^SN\t", keys[[nm]], ":\t"), sn, value = TRUE)
    if (length(hit)) out[[nm]] <- suppressWarnings(as.numeric(strsplit(hit[[1]], "\t", fixed = TRUE)[[1]][[3]]))
  }
  out
}

extract_bowtie2_log <- function(path) {
  lines <- read_lines_safe(path)
  if (length(lines) == 0) {
    return(c(aligner_input_reads = NA, unaligned_reads = NA, unique_aligned_reads = NA,
             multi_aligned_reads = NA, overall_alignment_pct = NA))
  }
  total_line <- grep(" reads; of these:$", lines, value = TRUE)
  total <- if (length(total_line)) suppressWarnings(as.numeric(sub("^\\s*([0-9,]+).*", "\\1", gsub(",", "", total_line[[1]])))) else NA_real_
  unaligned_line <- grep("aligned 0 times", lines, value = TRUE)
  unique_line <- grep("aligned exactly 1 time", lines, value = TRUE)
  multi_line <- grep("aligned >1 times", lines, value = TRUE)
  rate_line <- grep("overall alignment rate", lines, value = TRUE)
  parse_first_int <- function(x) if (length(x)) suppressWarnings(as.numeric(sub("^\\s*([0-9,]+).*", "\\1", gsub(",", "", x[[1]])))) else NA_real_
  rate <- if (length(rate_line)) suppressWarnings(as.numeric(sub("^\\s*([^%]+)%.*", "\\1", rate_line[[1]]))) else NA_real_
  c(
    aligner_input_reads = total,
    unaligned_reads = parse_first_int(unaligned_line),
    unique_aligned_reads = parse_first_int(unique_line),
    multi_aligned_reads = parse_first_int(multi_line),
    overall_alignment_pct = rate
  )
}

json_number <- function(text, path) {
  pattern <- paste0('"', paste(path, collapse = '"[[:space:]]*:[[:space:]]*\\{[^}]*"'), '"[[:space:]]*:[[:space:]]*([0-9.]+)')
  hit <- regexpr(pattern, text, perl = TRUE)
  if (hit[[1]] < 0) return(NA_real_)
  value <- regmatches(text, hit)
  suppressWarnings(as.numeric(sub(".*:[[:space:]]*([0-9.]+)$", "\\1", value)))
}

extract_fastp_json <- function(path) {
  existing <- first_existing(path)
  keys <- c(raw_reads = NA_real_, raw_bases = NA_real_, trimmed_reads = NA_real_, trimmed_bases = NA_real_, q30_rate_after = NA_real_)
  if (is.na(existing) || file.info(existing)$size == 0) return(keys)
  text <- paste(read_lines_safe(existing), collapse = "")
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    x <- tryCatch(jsonlite::fromJSON(existing), error = function(e) NULL)
    if (!is.null(x)) {
      keys[["raw_reads"]] <- num_or_na(x$summary$before_filtering$total_reads)
      keys[["raw_bases"]] <- num_or_na(x$summary$before_filtering$total_bases)
      keys[["trimmed_reads"]] <- num_or_na(x$summary$after_filtering$total_reads)
      keys[["trimmed_bases"]] <- num_or_na(x$summary$after_filtering$total_bases)
      keys[["q30_rate_after"]] <- num_or_na(x$summary$after_filtering$q30_rate) * 100
      return(keys)
    }
  }
  keys[["raw_reads"]] <- json_number(text, c("before_filtering", "total_reads"))
  keys[["raw_bases"]] <- json_number(text, c("before_filtering", "total_bases"))
  keys[["trimmed_reads"]] <- json_number(text, c("after_filtering", "total_reads"))
  keys[["trimmed_bases"]] <- json_number(text, c("after_filtering", "total_bases"))
  q30 <- json_number(text, c("after_filtering", "q30_rate"))
  keys[["q30_rate_after"]] <- ifelse(is.na(q30), NA_real_, q30 * 100)
  keys
}

find_peak_file <- function(peaks_dir, sample_id) {
  files <- list.files(file.path(peaks_dir, sample_id), pattern = "_peaks\\.(narrowPeak|broadPeak)$", full.names = TRUE)
  if (length(files) == 0) return(NA_character_)
  files[[1]]
}

peak_type_from_file <- function(path) {
  if (is.na(path)) return(NA_character_)
  if (grepl("\\.broadPeak$", path)) "broad" else "narrow"
}

markdown_table <- function(x, con, max_rows = Inf) {
  if (nrow(x) == 0) {
    writeLines("No rows available.", con)
    return(invisible(NULL))
  }
  y <- x
  if (is.finite(max_rows) && nrow(y) > max_rows) y <- head(y, max_rows)
  y[] <- lapply(y, function(col) {
    col <- as.character(col)
    col[is.na(col)] <- "NA"
    gsub("\\|", "/", col)
  })
  writeLines(paste0("| ", paste(names(y), collapse = " | "), " |"), con)
  writeLines(paste0("| ", paste(rep("---", ncol(y)), collapse = " | "), " |"), con)
  apply(y, 1, function(row) writeLines(paste0("| ", paste(row, collapse = " | "), " |"), con))
  if (is.finite(max_rows) && nrow(x) > max_rows) {
    writeLines(paste0("\nShowing first ", max_rows, " of ", nrow(x), " rows."), con)
  }
}

summarise_diff <- function(diff_results) {
  if (nrow(diff_results) == 0) {
    return(list(overall = data.frame(), by_contrast = data.frame(), by_method = data.frame()))
  }
  padj <- if ("padj" %in% names(diff_results)) suppressWarnings(as.numeric(diff_results$padj)) else rep(NA_real_, nrow(diff_results))
  pvalue <- if ("pvalue" %in% names(diff_results)) suppressWarnings(as.numeric(diff_results$pvalue)) else rep(NA_real_, nrow(diff_results))
  overall <- data.frame(
    tested_rows = nrow(diff_results),
    rows_with_pvalue = sum(!is.na(pvalue)),
    rows_with_padj = sum(!is.na(padj)),
    significant_padj_0_05 = sum(padj < 0.05, na.rm = TRUE),
    significant_padj_0_10 = sum(padj < 0.10, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  if ("contrast" %in% names(diff_results)) {
    by_contrast <- do.call(rbind, lapply(split(seq_len(nrow(diff_results)), diff_results$contrast), function(idx) {
      p <- padj[idx]
      data.frame(
        contrast = as.character(diff_results$contrast[idx[[1]]]),
        tested_rows = length(idx),
        rows_with_padj = sum(!is.na(p)),
        significant_padj_0_05 = sum(p < 0.05, na.rm = TRUE),
        significant_padj_0_10 = sum(p < 0.10, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    by_contrast <- by_contrast[order(by_contrast$significant_padj_0_05, decreasing = TRUE), , drop = FALSE]
  } else {
    by_contrast <- data.frame()
  }
  if ("method" %in% names(diff_results)) {
    by_method <- as.data.frame(table(diff_results$method), stringsAsFactors = FALSE)
    names(by_method) <- c("method", "rows")
  } else {
    by_method <- data.frame()
  }
  list(overall = overall, by_contrast = by_contrast, by_method = by_method)
}

args <- parse_args()
metadata <- read_tsv_safe(args[["metadata"]])
outdir <- args[["output-dir"]]
report <- args[["report"]]
report_dir <- dirname(report)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
table_suffix <- Sys.getenv("PIPELINE_TABLE_SUFFIX", unset = "")

required_cols <- c("sample_id", "condition", "mark_or_factor", "is_control")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) {
  stop("Metadata is missing required columns for report: ", paste(missing_cols, collapse = ", "))
}

is_control <- tolower(as.character(metadata$is_control)) %in% c("true", "1", "yes", "y")
sample_ids <- metadata$sample_id

fastp_metrics <- t(vapply(sample_ids, function(sid) {
  extract_fastp_json(file.path(outdir, "040-trimming", sid, paste0(sid, ".fastp.json")))
}, numeric(5)))

aligner_metrics <- t(vapply(sample_ids, function(sid) {
  logs <- list.files(file.path(outdir, "050-alignment", sid), pattern = "\\.(bowtie2|bwa)\\.log$", full.names = TRUE)
  if (length(logs) == 0) {
    return(c(aligner_input_reads = NA, unaligned_reads = NA, unique_aligned_reads = NA,
             multi_aligned_reads = NA, overall_alignment_pct = NA))
  }
  if (grepl("\\.bowtie2\\.log$", logs[[1]])) extract_bowtie2_log(logs[[1]]) else c(
    aligner_input_reads = NA, unaligned_reads = NA, unique_aligned_reads = NA,
    multi_aligned_reads = NA, overall_alignment_pct = NA
  )
}, numeric(5)))

final_flagstat <- t(vapply(sample_ids, function(sid) {
  extract_flagstat(file.path(outdir, "060-filtering", sid, paste0(sid, ".flagstat.txt")))
}, numeric(5)))

aligned_flagstat <- t(vapply(sample_ids, function(sid) {
  extract_flagstat(file.path(outdir, "050-alignment", sid, paste0(sid, ".aligned.flagstat.txt")))
}, numeric(5)))

final_stats <- t(vapply(sample_ids, function(sid) {
  extract_samtools_stats(file.path(outdir, "060-filtering", sid, paste0(sid, ".stats.txt")))
}, numeric(7)))

peak_files <- vapply(sample_ids, function(sid) find_peak_file(file.path(outdir, "080-peak-calling"), sid), character(1))
peak_counts <- vapply(peak_files, count_lines_safe, integer(1))
peak_types <- vapply(peak_files, peak_type_from_file, character(1))

sample_summary <- data.frame(
  sample_id = sample_ids,
  condition = metadata$condition,
  mark_or_factor = metadata$mark_or_factor,
  replicate = if ("replicate" %in% names(metadata)) metadata$replicate else NA,
  is_control = is_control,
  control_id = if ("control_id" %in% names(metadata)) metadata$control_id else NA,
  raw_reads_fastp = fastp_metrics[, "raw_reads"],
  trimmed_reads_fastp = fastp_metrics[, "trimmed_reads"],
  trim_retained_pct = ifelse(!is.na(fastp_metrics[, "raw_reads"]) & fastp_metrics[, "raw_reads"] > 0,
                             fastp_metrics[, "trimmed_reads"] / fastp_metrics[, "raw_reads"] * 100, NA_real_),
  q30_after_trim_pct = fastp_metrics[, "q30_rate_after"],
  aligner_input_reads = aligner_metrics[, "aligner_input_reads"],
  aligner_overall_alignment_pct = aligner_metrics[, "overall_alignment_pct"],
  aligned_bam_reads = aligned_flagstat[, "total"],
  aligned_bam_mapped_pct = aligned_flagstat[, "mapped_pct"],
  aligned_bam_duplicate_pct = aligned_flagstat[, "duplicate_pct"],
  final_filtered_reads = final_flagstat[, "total"],
  final_mapped_pct = final_flagstat[, "mapped_pct"],
  final_duplicate_pct = final_flagstat[, "duplicate_pct"],
  final_average_read_length = final_stats[, "average_length"],
  final_insert_size_average = final_stats[, "insert_size_average"],
  peak_type = peak_types,
  peak_count = peak_counts,
  peak_file = peak_files,
  filtered_bam = file.exists(file.path(outdir, "060-filtering", sample_ids, paste0(sample_ids, ".filtered.bam"))),
  stringsAsFactors = FALSE
)

annotation_summary <- read_tsv_safe(file.path(outdir, "090-peak-annotation", "annotation_summary.tsv"))
diff_results <- read_tsv_safe(file.path(outdir, "120-differential-binding", "differential_binding_results.tsv"))
diff_summary <- summarise_diff(diff_results)

consensus_files <- list.files(file.path(outdir, "110-consensus-peaks", "groups"), pattern = "\\.consensus\\.bed$", full.names = TRUE)
consensus_summary <- if (length(consensus_files) > 0) {
  data.frame(
    peak_set = sub("\\.consensus\\.bed$", "", basename(consensus_files)),
    consensus_peak_count = vapply(consensus_files, count_lines_safe, integer(1)),
    stringsAsFactors = FALSE
  )
} else {
  data.frame()
}

count_files <- list.files(file.path(outdir, "110-consensus-peaks", "counts"), pattern = "\\.counts\\.tsv(\\.gz)?$", full.names = TRUE)
count_summary <- if (length(count_files) > 0) {
  data.frame(
    count_matrix = basename(count_files),
    rows = vapply(count_files, count_lines_safe, integer(1)),
    stringsAsFactors = FALSE
  )
} else {
  data.frame()
}

annotation_totals <- if (nrow(annotation_summary) > 0 && all(c("class", "n") %in% names(annotation_summary))) {
  aggregate(n ~ class, annotation_summary, sum)
} else {
  data.frame()
}

group_summary <- aggregate(
  sample_id ~ condition + mark_or_factor,
  sample_summary[!sample_summary$is_control, , drop = FALSE],
  length
)
names(group_summary)[names(group_summary) == "sample_id"] <- "ip_replicates"
group_peaks <- aggregate(
  peak_count ~ condition + mark_or_factor,
  sample_summary[!sample_summary$is_control, , drop = FALSE],
  function(x) paste(c("min" = min(x, na.rm = TRUE), "median" = median(x, na.rm = TRUE), "max" = max(x, na.rm = TRUE)), collapse = " / ")
)
names(group_peaks)[names(group_peaks) == "peak_count"] <- "peak_count_min_median_max"
group_summary <- merge(group_summary, group_peaks, by = c("condition", "mark_or_factor"), all.x = TRUE)

recommendations <- character()
if (sum(sample_summary$is_control) == 0) {
  recommendations <- c(recommendations, "No input/control samples were provided. MACS peak calls should be interpreted as no-control calls with higher background risk.")
}
missing_control <- sample_summary$sample_id[!sample_summary$is_control & (is.na(sample_summary$control_id) | sample_summary$control_id == "")]
if (length(missing_control) > 0 && sum(sample_summary$is_control) > 0) {
  recommendations <- c(recommendations, paste("Some IP samples have no matched control_id:", paste(missing_control, collapse = ", ")))
}
low_mapping <- sample_summary$sample_id[!is.na(sample_summary$aligner_overall_alignment_pct) & sample_summary$aligner_overall_alignment_pct < 70]
if (length(low_mapping) > 0) {
  recommendations <- c(recommendations, paste("Inspect low aligner mapping rate samples:", paste(low_mapping, collapse = ", ")))
}
low_trim <- sample_summary$sample_id[!is.na(sample_summary$trim_retained_pct) & sample_summary$trim_retained_pct < 70]
if (length(low_trim) > 0) {
  recommendations <- c(recommendations, paste("Inspect samples with low read retention after trimming:", paste(low_trim, collapse = ", ")))
}
zero_peaks <- sample_summary$sample_id[!is.na(sample_summary$peak_count) & sample_summary$peak_count == 0 & !sample_summary$is_control]
if (length(zero_peaks) > 0) {
  recommendations <- c(recommendations, paste("No peaks detected for:", paste(zero_peaks, collapse = ", ")))
}
low_peaks <- sample_summary$sample_id[!is.na(sample_summary$peak_count) & sample_summary$peak_count > 0 & sample_summary$peak_count < 100 & !sample_summary$is_control]
if (length(low_peaks) > 0) {
  recommendations <- c(recommendations, paste("Very low peak counts (<100) detected:", paste(low_peaks, collapse = ", ")))
}
if (nrow(group_summary) > 0 && any(group_summary$ip_replicates < 2)) {
  low_rep <- paste(group_summary$condition[group_summary$ip_replicates < 2], group_summary$mark_or_factor[group_summary$ip_replicates < 2], sep = "/")
  recommendations <- c(recommendations, paste("Differential binding is underpowered for groups with <2 IP replicates:", paste(low_rep, collapse = ", ")))
}
if (nrow(diff_summary$overall) > 0 && diff_summary$overall$rows_with_padj[[1]] == 0) {
  recommendations <- c(recommendations, "Differential tables do not contain adjusted p-values. Check whether DESeq2 was available or whether fallback mode was used.")
}
if (length(recommendations) == 0) {
  recommendations <- "No automatic quality warnings were triggered."
}

sample_summary_file <- file.path(report_dir, paste0("chipseq_sample_qc_summary.tsv", table_suffix))
group_summary_file <- file.path(report_dir, paste0("chipseq_group_summary.tsv", table_suffix))
diff_summary_file <- file.path(report_dir, paste0("chipseq_differential_summary.tsv", table_suffix))
write_tsv_safe(sample_summary, sample_summary_file)
write_tsv_safe(group_summary, group_summary_file)
if (nrow(diff_summary$by_contrast) > 0) write_tsv_safe(diff_summary$by_contrast, diff_summary_file)

con <- file(report, open = "w", encoding = "UTF-8")
on.exit(close(con), add = TRUE)

writeLines(c(
  "# ChIP-seq final report",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Work directory:", outdir),
  "",
  "## Run overview",
  "",
  paste("Total samples:", nrow(metadata)),
  paste("IP samples:", sum(!is_control)),
  paste("Control/input samples:", sum(is_control)),
  paste("Conditions:", paste(sort(unique(metadata$condition)), collapse = ", ")),
  paste("Marks/factors:", paste(sort(unique(metadata$mark_or_factor)), collapse = ", ")),
  "",
  "## Important interpretation notes",
  "",
  "- `aligner_overall_alignment_pct` comes from the aligner log when available and is the preferred mapping-rate metric.",
  "- `aligned_bam_*` metrics come from the sorted BAM before MAPQ filtering, duplicate removal, and blacklist filtering when that file is available.",
  "- `final_mapped_pct` and `final_duplicate_pct` are calculated after filtering and duplicate removal; they describe the final BAM, not the original library complexity.",
  "- Differential binding counts below separate tested rows from statistically significant rows.",
  ""
), con)

writeLines("## Group summary\n", con)
markdown_table(group_summary, con, max_rows = 100)

writeLines(c("", "## Sample QC summary", ""), con)
sample_display <- sample_summary[, c("sample_id", "condition", "mark_or_factor", "is_control", "raw_reads_fastp",
                                     "trimmed_reads_fastp", "trim_retained_pct", "aligner_overall_alignment_pct",
                                     "aligned_bam_reads", "aligned_bam_mapped_pct", "final_filtered_reads",
                                     "peak_type", "peak_count"), drop = FALSE]
sample_display$trim_retained_pct <- percent(sample_display$trim_retained_pct)
sample_display$aligner_overall_alignment_pct <- percent(sample_display$aligner_overall_alignment_pct)
sample_display$aligned_bam_mapped_pct <- percent(sample_display$aligned_bam_mapped_pct)
markdown_table(sample_display, con, max_rows = 200)
writeLines(paste0("\nFull sample QC table: `", basename(sample_summary_file), "`"), con)

writeLines(c("", "## Peaks and consensus", ""), con)
writeLines(paste("Individual peak sets found:", sum(!is.na(sample_summary$peak_count))), con)
writeLines(paste("Total individual peaks:", fmt_num(sum(sample_summary$peak_count, na.rm = TRUE))), con)
if (nrow(consensus_summary) > 0) {
  writeLines(paste("Consensus peak sets:", nrow(consensus_summary)), con)
  markdown_table(consensus_summary[order(consensus_summary$consensus_peak_count, decreasing = TRUE), , drop = FALSE], con, max_rows = 30)
} else {
  writeLines("No consensus peak BED files were found.", con)
}
if (nrow(count_summary) > 0) {
  writeLines(c("", "Consensus count matrices:"), con)
  markdown_table(count_summary, con, max_rows = 50)
}

writeLines(c("", "## Peak annotation", ""), con)
if (nrow(annotation_totals) > 0) {
  markdown_table(annotation_totals[order(annotation_totals$n, decreasing = TRUE), , drop = FALSE], con)
  writeLines("", con)
}
if (nrow(annotation_summary) > 0) {
  writeLines("Annotation summary by peak set:", con)
  markdown_table(annotation_summary, con, max_rows = 80)
} else {
  writeLines("Peak annotation summary was not found.", con)
}

writeLines(c("", "## Differential binding", ""), con)
if (nrow(diff_summary$overall) > 0) {
  markdown_table(diff_summary$overall, con)
  if (nrow(diff_summary$by_method) > 0) {
    writeLines(c("", "Methods used:"), con)
    markdown_table(diff_summary$by_method, con)
  }
  if (nrow(diff_summary$by_contrast) > 0) {
    writeLines(c("", "Differential summary by contrast:"), con)
    markdown_table(diff_summary$by_contrast, con, max_rows = 50)
    writeLines(paste0("\nFull differential summary: `", basename(diff_summary_file), "`"), con)
  }
} else {
  writeLines("No differential binding results were generated.", con)
}

writeLines(c(
  "",
  "## Main output directories",
  "",
  "- `030-qc-fastq`: raw and post-trim FastQC/MultiQC",
  "- `040-trimming`: trimmed FASTQs and fastp reports",
  "- `050-alignment`: sorted BAMs and aligner logs",
  "- `060-filtering`: final filtered BAMs and BAM metrics",
  "- `070-qc-alignment`: alignment QC and MultiQC",
  "- `080-peak-calling`: MACS peak calls",
  "- `090-peak-annotation`: annotated peak tables",
  "- `100-tracks`: bigWig tracks",
  "- `110-consensus-peaks`: merged peaks and count matrices",
  "- `120-differential-binding`: differential enrichment tables and plots",
  "",
  "## Automatic recommendations",
  ""
), con)
writeLines(paste("-", recommendations), con)

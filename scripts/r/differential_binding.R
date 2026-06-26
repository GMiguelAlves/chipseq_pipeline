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

open_table <- function(path, mode = "rt") {
  if (grepl("\\.gz$", path)) gzfile(path, mode) else file(path, mode)
}

read_table_safe <- function(path, header = TRUE) {
  con <- open_table(path, "rt")
  on.exit(close(con), add = TRUE)
  read.table(con, sep = "\t", header = header, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
}

write_table_safe <- function(x, path, header = TRUE) {
  con <- open_table(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = header)
}

read_metadata <- function(path) {
  x <- read_table_safe(path, header = TRUE)
  x$is_control_norm <- tolower(as.character(x$is_control)) %in% c("true", "1", "yes", "y")
  x
}

contrast_pairs <- function(conditions, contrast_text) {
  conditions <- unique(as.character(conditions))
  if (!is.null(contrast_text) && nzchar(contrast_text)) {
    raw <- unlist(strsplit(contrast_text, ","))
    pairs <- strsplit(raw, ":")
    return(Filter(function(x) length(x) == 2, pairs))
  }
  if (length(conditions) < 2) return(list())
  combn(conditions, 2, simplify = FALSE)
}

plot_basic_qc <- function(counts, sample_info, out_prefix) {
  log_counts <- log2(counts + 1)
  if (ncol(log_counts) >= 2 && nrow(log_counts) >= 2) {
    pdf(paste0(out_prefix, ".pca.pdf"))
    pca <- prcomp(t(log_counts), scale. = TRUE)
    plot(pca$x[, 1], pca$x[, 2], pch = 19, col = as.integer(factor(sample_info$condition)),
         xlab = "PC1", ylab = "PC2", main = "Peak count PCA")
    text(pca$x[, 1], pca$x[, 2], labels = sample_info$sample_id, pos = 3, cex = 0.7)
    legend("topright", legend = levels(factor(sample_info$condition)), col = seq_along(levels(factor(sample_info$condition))), pch = 19, cex = 0.7)
    dev.off()

    top_var <- order(apply(log_counts, 1, var), decreasing = TRUE)
    top_var <- head(top_var, min(500, length(top_var)))
    pdf(paste0(out_prefix, ".heatmap.pdf"))
    heatmap(as.matrix(log_counts[top_var, , drop = FALSE]), Colv = NA, scale = "row", margins = c(8, 5))
    dev.off()
  }
}

fallback_contrast <- function(counts, sample_info, a, b) {
  group_a <- sample_info$sample_id[sample_info$condition == a]
  group_b <- sample_info$sample_id[sample_info$condition == b]
  mean_a <- rowMeans(log2(counts[, group_a, drop = FALSE] + 1))
  mean_b <- rowMeans(log2(counts[, group_b, drop = FALSE] + 1))
  data.frame(
    peak_id = rownames(counts),
    baseMean = rowMeans(counts),
    log2FC = mean_a - mean_b,
    pvalue = NA_real_,
    padj = NA_real_,
    contrast = paste(a, "vs", b, sep = "_"),
    method = "fallback_log2_mean_difference",
    stringsAsFactors = FALSE
  )
}

run_deseq2 <- function(counts, sample_info, pairs) {
  suppressPackageStartupMessages(library(DESeq2))
  condition <- factor(sample_info$condition)
  coldata <- data.frame(row.names = sample_info$sample_id, condition = condition)
  dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = coldata, design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  out <- list()
  for (pair in pairs) {
    res <- results(dds, contrast = c("condition", pair[[1]], pair[[2]]))
    tab <- as.data.frame(res)
    tab$peak_id <- rownames(tab)
    tab$contrast <- paste(pair[[1]], "vs", pair[[2]], sep = "_")
    tab$method <- "DESeq2"
    names(tab)[names(tab) == "log2FoldChange"] <- "log2FC"
    out[[length(out) + 1]] <- tab[, intersect(c("peak_id", "baseMean", "log2FC", "pvalue", "padj", "contrast", "method"), names(tab))]
  }
  do.call(rbind, out)
}

plot_contrast <- function(res, out_prefix) {
  if (nrow(res) == 0) return(invisible(NULL))
  pdf(paste0(out_prefix, ".ma.pdf"))
  plot(log10(res$baseMean + 1), res$log2FC, pch = 16, cex = 0.4, xlab = "log10 baseMean", ylab = "log2FC", main = "MA plot")
  abline(h = 0, col = "red")
  dev.off()

  if (all(is.na(res$pvalue))) return(invisible(NULL))
  pdf(paste0(out_prefix, ".volcano.pdf"))
  plot(res$log2FC, -log10(res$pvalue), pch = 16, cex = 0.4, xlab = "log2FC", ylab = "-log10 pvalue", main = "Volcano plot")
  abline(v = c(-1, 1), col = "grey")
  dev.off()
}

args <- parse_args()
required <- c("metadata", "counts-dir", "outdir", "min-replicates")
missing <- required[!vapply(required, function(x) !is.null(args[[x]]) && nzchar(args[[x]]), logical(1))]
if (length(missing) > 0) stop("Missing required argument(s): ", paste(missing, collapse = ", "))

dir.create(args[["outdir"]], recursive = TRUE, showWarnings = FALSE)
table_suffix <- if (!is.null(args[["table-suffix"]])) args[["table-suffix"]] else Sys.getenv("PIPELINE_TABLE_SUFFIX", unset = "")
metadata <- read_metadata(args[["metadata"]])
ip_metadata <- metadata[!metadata$is_control_norm, , drop = FALSE]
sample_order <- ip_metadata$sample_id
min_reps <- as.integer(args[["min-replicates"]])

count_files <- list.files(args[["counts-dir"]], pattern = "\\.counts\\.tsv(\\.gz)?$", full.names = TRUE)
if (length(count_files) == 0) {
  writeLines("No consensus count files found.", file.path(args[["outdir"]], "differential_binding_skipped.txt"))
  quit(save = "no")
}

all_results <- list()
for (count_file in count_files) {
  raw <- read_table_safe(count_file, header = FALSE)
  if (ncol(raw) < 4) next
  n_counts <- ncol(raw) - 3
  samples <- head(sample_order, n_counts)
  counts <- as.matrix(raw[, 4:ncol(raw), drop = FALSE])
  storage.mode(counts) <- "numeric"
  colnames(counts) <- samples
  rownames(counts) <- paste(raw[[1]], raw[[2]], raw[[3]], sep = ":")
  sample_info <- ip_metadata[match(samples, ip_metadata$sample_id), , drop = FALSE]
  sample_info <- sample_info[!is.na(sample_info$sample_id), , drop = FALSE]
  counts <- counts[, sample_info$sample_id, drop = FALSE]

  reps <- table(sample_info$condition)
  eligible_conditions <- names(reps)[reps >= min_reps]
  sample_info <- sample_info[sample_info$condition %in% eligible_conditions, , drop = FALSE]
  counts <- counts[, sample_info$sample_id, drop = FALSE]
  if (length(unique(sample_info$condition)) < 2) {
    next
  }

  base_name <- sub("\\.gz$", "", basename(count_file))
  prefix <- file.path(args[["outdir"]], sub("\\.counts\\.tsv$", "", base_name))
  plot_basic_qc(counts, sample_info, prefix)
  pairs <- contrast_pairs(sample_info$condition, args[["contrasts"]])
  pairs <- Filter(function(p) all(p %in% sample_info$condition), pairs)
  if (length(pairs) == 0) next

  if (requireNamespace("DESeq2", quietly = TRUE)) {
    res <- run_deseq2(counts, sample_info, pairs)
  } else {
    res <- do.call(rbind, lapply(pairs, function(p) fallback_contrast(counts, sample_info, p[[1]], p[[2]])))
  }
  write_table_safe(res, paste0(prefix, ".differential.tsv", table_suffix))
  split_res <- split(res, res$contrast)
  for (nm in names(split_res)) {
    plot_contrast(split_res[[nm]], paste0(prefix, ".", nm))
  }
  all_results[[length(all_results) + 1]] <- res
}

combined <- if (length(all_results) > 0) do.call(rbind, all_results) else data.frame()
write_table_safe(combined, file.path(args[["outdir"]], paste0("differential_binding_results.tsv", table_suffix)))
if (nrow(combined) == 0) {
  writeLines("No eligible differential contrasts were found. Check conditions, contrasts, and replicate counts.", file.path(args[["outdir"]], "differential_binding_skipped.txt"))
}

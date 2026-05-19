#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ShortRead)
})

source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]])), "lib_pipeline.R"))

args <- parse_cli(list(run_id = "run_01_default", reads_per_file = "1000"))
project <- script_project_root()
run_id <- args$run_id
reads_per_file <- as.integer(args$reads_per_file)

report_dir <- file.path(project, "reports", "preprocessing_summary")
qc_dir <- file.path(project, "results", "08_qc_tables")
ensure_dir(report_dir, qc_dir)

metadata <- read_run_metadata(project, run_id)
cutadapt_stats <- read_tsv(file.path(run_dir(project, "01_cutadapt", run_id), "cutadapt_stats.tsv"))
filter_stats <- read_tsv(file.path(run_dir(project, "02_filter_trim", run_id), "dada2_filter_stats.tsv"))
merge_stats <- read_tsv(file.path(run_dir(project, "05_chimera", run_id), "dada2_merge_chimera_stats.tsv"))
seqtab <- readRDS(file.path(run_dir(project, "04_denoise_merge", run_id), "seqtab.rds"))
seqtab_nochim <- readRDS(file.path(run_dir(project, "05_chimera", run_id), "seqtab_nochim.rds"))

track <- merge(metadata[, .(sample_id, sample_type, include_biological_analysis, diet_code, timepoint_day, group)],
               cutadapt_stats, by = "sample_id", all.x = TRUE)
track <- merge(track, filter_stats[, .(sample_id, reads.in, reads.out)], by = "sample_id", all.x = TRUE)
track <- merge(track, merge_stats, by = "sample_id", all.x = TRUE, suffixes = c("", "_dada2"))
track[, `:=`(
  cutadapt_retention_pct = pairs_written / pairs_processed * 100,
  raw_to_filtered_pct = reads.out / pairs_processed * 100,
  raw_to_nonchim_pct = nonchim / pairs_processed * 100,
  denoise_forward_loss_pct = (reads.out - denoised_forward) / reads.out * 100,
  denoise_reverse_loss_pct = (reads.out - denoised_reverse) / reads.out * 100,
  merge_from_filtered_output_pct = merged / reads.out * 100,
  merge_loss_from_denoised_pct = (pmin(denoised_forward, denoised_reverse) - merged) / pmin(denoised_forward, denoised_reverse) * 100,
  chimera_loss_of_merged_pct = (merged - nonchim) / merged * 100
)]

sample_first_reads <- function(file, n) {
  if (!file.exists(file)) return(NULL)
  stream <- FastqStreamer(file, n = n)
  on.exit(close(stream))
  fq <- yield(stream)
  if (length(fq) == 0) return(NULL)
  q <- as(quality(fq), "matrix")
  len <- width(sread(fq))
  list(q = q, length = len)
}

quality_by_position <- function(files, label, n) {
  pos_tables <- list()
  length_rows <- list()
  for (i in seq_along(files)) {
    sampled <- sample_first_reads(files[[i]], n)
    if (is.null(sampled)) next
    q <- sampled$q
    length_rows[[length(length_rows) + 1L]] <- data.table(
      file = basename(files[[i]]),
      reads_sampled = nrow(q),
      min_length = min(sampled$length),
      median_length = median(sampled$length),
      max_length = max(sampled$length)
    )
    pos_tables[[length(pos_tables) + 1L]] <- data.table(
      position = seq_len(ncol(q)),
      mean_q = colMeans(q),
      median_q = apply(q, 2, median),
      q10 = apply(q, 2, quantile, probs = 0.10, na.rm = TRUE),
      q25 = apply(q, 2, quantile, probs = 0.25, na.rm = TRUE),
      reads_at_position = colSums(!is.na(q))
    )
  }
  if (length(pos_tables) == 0) {
    return(list(quality = data.table(), lengths = data.table()))
  }
  qdt <- rbindlist(pos_tables)
  qdt <- qdt[, .(
    mean_q = mean(mean_q, na.rm = TRUE),
    median_q = median(median_q, na.rm = TRUE),
    q10 = median(q10, na.rm = TRUE),
    q25 = median(q25, na.rm = TRUE),
    reads_at_position = sum(reads_at_position, na.rm = TRUE)
  ), by = position]
  qdt[, dataset := label]
  list(quality = qdt, lengths = rbindlist(length_rows, fill = TRUE)[, dataset := label])
}

log_message("Sampling qualities for", run_id, "with", reads_per_file, "reads/file")
raw_r1 <- quality_by_position(metadata$forward_symlink, "raw_R1", reads_per_file)
raw_r2 <- quality_by_position(metadata$reverse_symlink, "raw_R2", reads_per_file)
trim_r1 <- quality_by_position(cutadapt_stats$r1_trimmed, "trimmed_R1", reads_per_file)
trim_r2 <- quality_by_position(cutadapt_stats$r2_trimmed, "trimmed_R2", reads_per_file)
qual <- rbindlist(list(raw_r1$quality, raw_r2$quality, trim_r1$quality, trim_r2$quality), fill = TRUE)
len_samples <- rbindlist(list(raw_r1$lengths, raw_r2$lengths, trim_r1$lengths, trim_r2$lengths), fill = TRUE)

write_tsv(track, file.path(qc_dir, paste0(run_id, "_optimization_track_reads.tsv")))
write_tsv(qual, file.path(qc_dir, paste0(run_id, "_quality_by_position_sampled.tsv")))
write_tsv(len_samples, file.path(qc_dir, paste0(run_id, "_read_length_sampled.tsv")))

qual_plot <- ggplot(qual, aes(position, q10, color = dataset)) +
  geom_line(linewidth = 0.35) +
  geom_hline(yintercept = c(20, 25, 30), linetype = "dashed", color = "grey65") +
  labs(x = "Position", y = "Median per-file Q10", title = paste(run_id, "sampled quality tails")) +
  theme_minimal(base_size = 11)
ggsave(file.path(report_dir, paste0(run_id, "_quality_q10_by_position.pdf")), qual_plot, width = 10, height = 6)

asv_lengths_all <- data.table(sequence = colnames(seqtab), length = nchar(colnames(seqtab)), abundance = colSums(seqtab), nonchim = FALSE)
asv_lengths_nochim <- data.table(sequence = colnames(seqtab_nochim), length = nchar(colnames(seqtab_nochim)), abundance = colSums(seqtab_nochim), nonchim = TRUE)
asv_len_summary <- rbindlist(list(
  asv_lengths_all[, .(stage = "merged", asvs = .N, reads = sum(abundance), min = min(length), q25 = quantile(length, .25), median = median(length), q75 = quantile(length, .75), max = max(length))],
  asv_lengths_nochim[, .(stage = "nonchim", asvs = .N, reads = sum(abundance), min = min(length), q25 = quantile(length, .25), median = median(length), q75 = quantile(length, .75), max = max(length))]
))
write_tsv(asv_len_summary, file.path(qc_dir, paste0(run_id, "_asv_length_summary.tsv")))

summ <- function(x) {
  sprintf("min %.2f | q25 %.2f | median %.2f | mean %.2f | q75 %.2f | max %.2f",
          min(x, na.rm = TRUE), quantile(x, 0.25, na.rm = TRUE), median(x, na.rm = TRUE),
          mean(x, na.rm = TRUE), quantile(x, 0.75, na.rm = TRUE), max(x, na.rm = TRUE))
}

first_sustained_below <- function(position, value, threshold, window = 10L) {
  o <- order(position)
  position <- position[o]
  value <- value[o]
  below <- value < threshold
  if (length(below) < window) return(NA_integer_)
  runs <- stats::filter(as.integer(below), rep(1L, window), sides = 1)
  hit <- which(runs >= window)[1]
  if (is.na(hit)) NA_integer_ else position[[hit - window + 1L]]
}

q_tail <- qual[, .(
  first_pos_q10_below_30 = suppressWarnings(min(position[q10 < 30], na.rm = TRUE)),
  first_pos_q10_below_25 = suppressWarnings(min(position[q10 < 25], na.rm = TRUE)),
  first_pos_q10_below_20 = suppressWarnings(min(position[q10 < 20], na.rm = TRUE)),
  first_sustained_10bp_q10_below_30 = first_sustained_below(position, q10, 30),
  first_sustained_10bp_q10_below_25 = first_sustained_below(position, q10, 25),
  first_sustained_10bp_q10_below_20 = first_sustained_below(position, q10, 20),
  median_q10_from_200_to_end = median(q10[position >= 200], na.rm = TRUE),
  min_q10_from_200_to_end = min(q10[position >= 200], na.rm = TRUE),
  last_position = max(position)
), by = dataset]
for (col in c("first_pos_q10_below_30", "first_pos_q10_below_25", "first_pos_q10_below_20",
              "first_sustained_10bp_q10_below_30", "first_sustained_10bp_q10_below_25",
              "first_sustained_10bp_q10_below_20")) {
  q_tail[is.infinite(get(col)), (col) := NA_integer_]
}
write_tsv(q_tail, file.path(qc_dir, paste0(run_id, "_quality_threshold_positions.tsv")))

low_depth <- track[include_biological_analysis == "yes"][order(nonchim)][1:min(.N, 15)]
high_chim <- track[order(-chimera_loss_of_merged_pct)][1:min(.N, 15)]
low_merge <- track[order(merge_from_filtered_output_pct)][1:min(.N, 15)]
low_filter <- track[order(filter_retention_pct)][1:min(.N, 15)]

decision <- "Keep run_01_default as a strong baseline; test one adjusted DADA2 run aimed at recovering filter-stage losses without increasing chimera/merge penalties."
if (median(track$merge_from_filtered_output_pct, na.rm = TRUE) < 95) {
  decision <- "Prioritize truncLen tuning before relaxing maxEE, because merge retention suggests tail noise or insufficient overlap."
} else if (median(track$filter_retention_pct, na.rm = TRUE) < 88 && median(track$merge_from_filtered_output_pct, na.rm = TRUE) >= 95) {
  decision <- "Prioritize a modest maxEE relaxation; merge is robust enough that the main recoverable loss is filterAndTrim."
}

report <- c(
  paste0("# ", run_id, " DADA2 optimization evaluation"),
  "",
  "## Decision",
  "",
  paste0("- ", decision),
  "- Suggested next run: `run_02_dada2_maxee_3_3` with the same cutadapt outputs, `truncLen=c(0,0)`, `maxEE=c(3,3)`, `maxN=0`, `truncQ=2`, `mergePairs(maxMismatch=0)`, and the same chimera method. Compare against run_01 on raw-to-nonchim retention, chimera read loss, ASV length distribution, and taxonomy completeness.",
  "- Do not loosen `maxMismatch` yet: current merging is already high, so allowing mismatches would mostly reduce confidence rather than solve the main bottleneck.",
  "",
  "## Read Retention",
  "",
  paste0("- Cutadapt retention (%): ", summ(track$cutadapt_retention_pct)),
  paste0("- DADA2 filter retention from primer-trimmed reads (%): ", summ(track$filter_retention_pct)),
  paste0("- Merged reads relative to filterAndTrim output (%): ", summ(track$merge_from_filtered_output_pct)),
  paste0("- Merged reads relative to primer-trimmed input (%): ", summ(track$merge_retention_pct)),
  paste0("- Non-chimeric reads relative to primer-trimmed input (%): ", summ(track$nonchim_retention_pct)),
  paste0("- Raw-to-nonchim retention (%): ", summ(track$raw_to_nonchim_pct)),
  paste0("- Chimera read loss (% of merged reads): ", summ(track$chimera_loss_of_merged_pct)),
  "",
  "## Quality And Length",
  "",
  paste0("- Sampled reads per FASTQ for quality table: ", reads_per_file),
  paste0("- ASV length summary: ", paste(sprintf("%s median=%s q25=%s q75=%s min=%s max=%s ASVs=%s",
                                                asv_len_summary$stage, asv_len_summary$median,
                                                asv_len_summary$q25, asv_len_summary$q75,
                                                asv_len_summary$min, asv_len_summary$max,
                                                asv_len_summary$asvs), collapse = "; ")),
  "",
  "## Quality Threshold Positions",
  "",
  paste(capture.output(print(q_tail)), collapse = "\n"),
  "",
  "## Samples To Inspect",
  "",
  "### Lowest Final Depth",
  paste(capture.output(print(low_depth[, .(sample_id, group, pairs_processed, reads.out, merged, nonchim, raw_to_nonchim_pct)])), collapse = "\n"),
  "",
  "### Highest Chimera Loss",
  paste(capture.output(print(high_chim[, .(sample_id, group, merged, nonchim, chimera_loss_of_merged_pct)])), collapse = "\n"),
  "",
  "### Lowest Merge Retention",
  paste(capture.output(print(low_merge[, .(sample_id, group, reads.out, merged, merge_from_filtered_output_pct)])), collapse = "\n"),
  "",
  "### Lowest Filter Retention",
  paste(capture.output(print(low_filter[, .(sample_id, group, pairs_written, reads.out, filter_retention_pct)])), collapse = "\n"),
  "",
  "## Files",
  "",
  paste0("- `results/08_qc_tables/", run_id, "_optimization_track_reads.tsv`"),
  paste0("- `results/08_qc_tables/", run_id, "_quality_by_position_sampled.tsv`"),
  paste0("- `results/08_qc_tables/", run_id, "_quality_threshold_positions.tsv`"),
  paste0("- `results/08_qc_tables/", run_id, "_asv_length_summary.tsv`"),
  paste0("- `reports/preprocessing_summary/", run_id, "_quality_q10_by_position.pdf`")
)
writeLines(report, file.path(report_dir, paste0(run_id, "_optimization_evaluation.md")))
log_message("Optimization evaluation complete:", file.path(report_dir, paste0(run_id, "_optimization_evaluation.md")))

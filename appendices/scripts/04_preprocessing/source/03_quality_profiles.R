#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dada2)
})

source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]])), "lib_pipeline.R"))

args <- parse_cli(list(run_id = "run_01_default", n = "12"))
project <- script_project_root()
run_id <- args$run_id
n <- as.integer(args$n)
metadata <- read_run_metadata(project, run_id)
cut_stats <- read_tsv(file.path(run_dir(project, "01_cutadapt", run_id), "cutadapt_stats.tsv"))
ensure_dir(file.path(project, "reports", "preprocessing_summary"))

plot_pdf <- function(files, out, title) {
  files <- files[file.exists(files)]
  files <- files[seq_len(min(length(files), n))]
  if (length(files) == 0) return(invisible(FALSE))
  pdf(out, width = 12, height = 7)
  print(plotQualityProfile(files) + ggplot2::ggtitle(title))
  dev.off()
  invisible(TRUE)
}

plot_pdf(metadata$forward_symlink, file.path(project, "reports", "preprocessing_summary", paste0(run_id, "_raw_R1_quality.pdf")), "Raw R1 quality")
plot_pdf(metadata$reverse_symlink, file.path(project, "reports", "preprocessing_summary", paste0(run_id, "_raw_R2_quality.pdf")), "Raw R2 quality")
plot_pdf(cut_stats$r1_trimmed, file.path(project, "reports", "preprocessing_summary", paste0(run_id, "_trimmed_R1_quality.pdf")), "Cutadapt-trimmed R1 quality")
plot_pdf(cut_stats$r2_trimmed, file.path(project, "reports", "preprocessing_summary", paste0(run_id, "_trimmed_R2_quality.pdf")), "Cutadapt-trimmed R2 quality")
log_message("Quality profile PDFs complete for", run_id)

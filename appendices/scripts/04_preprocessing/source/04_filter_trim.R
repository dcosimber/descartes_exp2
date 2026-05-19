#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dada2)
})

source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]])), "lib_pipeline.R"))

args <- parse_cli(list(run_id = "run_01_default"))
project <- script_project_root()
run_id <- args$run_id
params <- load_run_params(project, run_id)
filt <- params$dada2$filter_and_trim

cut_stats <- read_tsv(file.path(run_dir(project, "01_cutadapt", run_id), "cutadapt_stats.tsv"))
outdir <- run_dir(project, "02_filter_trim", run_id)
fastq_dir <- file.path(outdir, "filtered_fastq")
ensure_dir(outdir, fastq_dir, file.path(project, "logs", "commands"), file.path(project, "logs", "stdout_stderr"))

filt_f <- file.path(fastq_dir, paste0(cut_stats$sample_id, "_R1.filtered.fastq.gz"))
filt_r <- file.path(fastq_dir, paste0(cut_stats$sample_id, "_R2.filtered.fastq.gz"))

cmd_log <- file.path(project, "logs", "commands", paste0("04_filter_trim_", run_id, ".commands.R"))
writeLines(capture.output(str(filt)), cmd_log)

log_message("Running DADA2 filterAndTrim for", nrow(cut_stats), "samples.")
filter_stats <- filterAndTrim(
  fwd = cut_stats$r1_trimmed,
  filt = filt_f,
  rev = cut_stats$r2_trimmed,
  filt.rev = filt_r,
  truncLen = c(as.integer(filt$trunc_len_f), as.integer(filt$trunc_len_r)),
  trimLeft = c(as.integer(filt$trim_left_f), as.integer(filt$trim_left_r)),
  maxN = as.integer(filt$max_n),
  maxEE = c(as.numeric(filt$max_ee_f), as.numeric(filt$max_ee_r)),
  truncQ = as.integer(filt$trunc_q),
  rm.phix = isTRUE(filt$rm_phix),
  compress = isTRUE(filt$compress),
  multithread = TRUE
)

filter_dt <- as.data.table(filter_stats, keep.rownames = "input_path")
filter_dt[, sample_id := cut_stats$sample_id]
filter_dt[, `:=`(filt_f = filt_f, filt_r = filt_r)]
setcolorder(filter_dt, c("sample_id", setdiff(names(filter_dt), "sample_id")))

write_tsv(filter_dt, file.path(outdir, "dada2_filter_stats.tsv"))
write_tsv(filter_dt, file.path(project, "results", "08_qc_tables", paste0(run_id, "_dada2_filter_stats.tsv")))
saveRDS(list(params = filt, stats = filter_dt), file.path(outdir, "filter_trim_params_and_stats.rds"))
log_message("filterAndTrim complete.")

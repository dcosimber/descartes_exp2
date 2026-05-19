#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]])), "lib_pipeline.R"))

args <- parse_cli(list(run_id = "run_01_default", sample_subset = ""))
project <- script_project_root()
run_id <- args$run_id
params <- load_run_params(project, run_id)
metadata <- read_run_metadata(project, run_id)
subset_ids <- split_subset(args$sample_subset)
if (length(subset_ids) > 0) metadata <- metadata[sample_id %in% subset_ids]

outdir <- run_dir(project, "01_cutadapt", run_id)
fastq_dir <- file.path(outdir, "trimmed_fastq")
sample_log_dir <- file.path(outdir, "sample_logs")
ensure_dir(outdir, fastq_dir, sample_log_dir, file.path(project, "logs", "commands"), file.path(project, "logs", "stdout_stderr"))

cmd_log <- file.path(project, "logs", "commands", paste0("02_run_cutadapt_", run_id, ".commands.sh"))
stdout_log <- file.path(project, "logs", "stdout_stderr", paste0("02_run_cutadapt_", run_id, ".log"))
if (file.exists(cmd_log)) unlink(cmd_log)
if (file.exists(stdout_log)) unlink(stdout_log)

cut <- params$cutadapt
primers <- params$primers
cutadapt_bin <- Sys.which("cutadapt")
if (cutadapt_bin == "") stop("cutadapt not found in PATH. Activate the conda environment first.")

parse_cutadapt_log <- function(file) {
  txt <- paste(readLines(file, warn = FALSE), collapse = "\n")
  grab <- function(pattern) {
    m <- regexpr(pattern, txt, perl = TRUE)
    if (m < 0) return(NA_real_)
    as.numeric(gsub(",", "", regmatches(txt, m)))
  }
  grab_pct <- function(pattern) {
    m <- regexec(pattern, txt, perl = TRUE)
    hit <- regmatches(txt, m)[[1]]
    if (length(hit) < 2) NA_real_ else as.numeric(hit[[2]])
  }
  data.table(
    pairs_processed = grab("(?<=Total read pairs processed:)\\s*[0-9,]+"),
    r1_with_adapter_pct = grab_pct("Read 1 with adapter:\\s*[0-9,]+ \\(([0-9.]+)%\\)"),
    r2_with_adapter_pct = grab_pct("Read 2 with adapter:\\s*[0-9,]+ \\(([0-9.]+)%\\)"),
    pairs_too_short_pct = grab_pct("Pairs that were too short:\\s*[0-9,]+ \\(([0-9.]+)%\\)"),
    pairs_discarded_untrimmed_pct = grab_pct("Pairs discarded as untrimmed:\\s*[0-9,]+ \\(([0-9.]+)%\\)"),
    pairs_written = grab("(?<=Pairs written \\(passing filters\\):)\\s*[0-9,]+"),
    pairs_written_pct = grab_pct("Pairs written \\(passing filters\\):\\s*[0-9,]+ \\(([0-9.]+)%\\)")
  )
}

stats <- list()
for (i in seq_len(nrow(metadata))) {
  sid <- metadata$sample_id[i]
  out_f <- file.path(fastq_dir, paste0(sid, "_R1.trimmed.fastq.gz"))
  out_r <- file.path(fastq_dir, paste0(sid, "_R2.trimmed.fastq.gz"))
  sample_log <- file.path(sample_log_dir, paste0(sid, ".cutadapt.log"))
  cmd <- c(
    "--cores", as.character(cut$cores),
    "--error-rate", as.character(cut$error_rate),
    "--times", as.character(cut$times),
    "--overlap", as.character(cut$overlap),
    "--minimum-length", as.character(cut$minimum_length),
    "-q", paste0(cut$quality_cutoff_forward, ",", cut$quality_cutoff_reverse),
    "-g", primers$forward,
    "-G", primers$reverse
  )
  if (isTRUE(cut$match_read_wildcards)) cmd <- c(cmd, "--match-read-wildcards")
  if (isTRUE(cut$discard_untrimmed)) cmd <- c(cmd, "--discard-untrimmed")
  cmd <- c(cmd, "-o", out_f, "-p", out_r, metadata$forward_symlink[i], metadata$reverse_symlink[i])
  write(paste(shQuote(cutadapt_bin), paste(shQuote(cmd), collapse = " ")), cmd_log, append = TRUE)
  log_message("Cutadapt", sid)
  rc <- system2(cutadapt_bin, cmd, stdout = sample_log, stderr = sample_log)
  if (!identical(rc, 0L)) stop("cutadapt failed for sample ", sid, ". See ", sample_log)
  s <- parse_cutadapt_log(sample_log)
  if (!isTRUE(cut$discard_untrimmed) && is.na(s$pairs_discarded_untrimmed_pct)) {
    s$pairs_discarded_untrimmed_pct <- 0
  }
  s[, `:=`(
    sample_id = sid,
    r1_trimmed = out_f,
    r2_trimmed = out_r,
    sample_log = sample_log
  )]
  stats[[sid]] <- s
}

stats_dt <- rbindlist(stats, fill = TRUE)
setcolorder(stats_dt, c("sample_id", setdiff(names(stats_dt), "sample_id")))
write_tsv(stats_dt, file.path(outdir, "cutadapt_stats.tsv"))
write_tsv(stats_dt, file.path(project, "results", "08_qc_tables", paste0(run_id, "_cutadapt_stats.tsv")))
log_message("Cutadapt complete for", nrow(stats_dt), "samples.")

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
learn <- params$dada2$learn_errors

filter_stats <- read_tsv(file.path(run_dir(project, "02_filter_trim", run_id), "dada2_filter_stats.tsv"))
outdir <- run_dir(project, "03_error_models", run_id)
ensure_dir(outdir, file.path(project, "reports", "preprocessing_summary"), file.path(project, "logs", "commands"))

cmd_log <- file.path(project, "logs", "commands", paste0("05_learn_errors_", run_id, ".commands.R"))
writeLines(capture.output(str(learn)), cmd_log)

log_message("Learning forward errors.")
err_f <- learnErrors(
  filter_stats$filt_f,
  multithread = TRUE,
  randomize = isTRUE(learn$randomize),
  nbases = as.numeric(learn$nbases)
)
log_message("Learning reverse errors.")
err_r <- learnErrors(
  filter_stats$filt_r,
  multithread = TRUE,
  randomize = isTRUE(learn$randomize),
  nbases = as.numeric(learn$nbases)
)

saveRDS(err_f, file.path(outdir, "errF.rds"))
saveRDS(err_r, file.path(outdir, "errR.rds"))

pdf(file.path(project, "reports", "preprocessing_summary", paste0(run_id, "_errF.pdf")), width = 10, height = 7)
print(plotErrors(err_f, nominalQ = TRUE))
dev.off()
pdf(file.path(project, "reports", "preprocessing_summary", paste0(run_id, "_errR.pdf")), width = 10, height = 7)
print(plotErrors(err_r, nominalQ = TRUE))
dev.off()

log_message("Error models complete.")

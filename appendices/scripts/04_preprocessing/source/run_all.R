#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]])), "lib_pipeline.R"))

args <- parse_cli(list(run_id = "run_01_default", sample_subset = ""))
project <- script_project_root()
run_id <- args$run_id
subset_arg <- args$sample_subset

scripts <- c(
  "01_build_metadata.R",
  "02_run_cutadapt.R",
  "03_quality_profiles.R",
  "04_filter_trim.R",
  "05_learn_errors.R",
  "06_denoise_merge_chimera.R",
  "07_assign_taxonomy.R",
  "08_export_outputs.R",
  "09_evaluate_run.R",
  "10_prepare_downstream_inputs.R"
)

ensure_dir(file.path(project, "logs", "stdout_stderr"))
log_file <- file.path(project, "logs", "stdout_stderr", paste0("run_all_", run_id, ".log"))
log_message("Starting full pipeline", run_id)

for (s in scripts) {
  script_path <- file.path(project, "scripts", s)
  cmd <- c(script_path, "--run-id", run_id)
  if (!is.null(subset_arg) && !isTRUE(subset_arg) && !is.na(subset_arg) && nzchar(subset_arg)) {
    cmd <- c(cmd, "--sample-subset", subset_arg)
  }
  log_message("Running", s)
  rc <- system2("Rscript", cmd, stdout = log_file, stderr = log_file)
  if (!identical(rc, 0L)) stop("Pipeline step failed: ", s, ". See ", log_file)
}

log_message("Full pipeline complete:", run_id)

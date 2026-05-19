#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(tools)
})

source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]])), "lib_pipeline.R"))

args <- parse_cli(list(run_id = "run_01_default", sample_subset = ""))
project <- script_project_root()
params <- load_run_params(project, args$run_id)

raw_dir <- normalizePath(file.path(project, params$project$raw_data_dir), mustWork = TRUE)
md5_file <- normalizePath(file.path(project, params$project$raw_md5_file), mustWork = TRUE)
metadata_dir <- file.path(project, "metadata")
symlink_root <- file.path(project, "data", "raw_reads", "by_sample")
ensure_dir(metadata_dir, symlink_root)

subset_ids <- split_subset(args$sample_subset)
is_subset_run <- length(subset_ids) > 0
log_message("Building metadata from", raw_dir)

r1 <- list.files(raw_dir, pattern = "\\.raw_1\\.fastq\\.gz$", recursive = TRUE, full.names = TRUE)
r2 <- list.files(raw_dir, pattern = "\\.raw_2\\.fastq\\.gz$", recursive = TRUE, full.names = TRUE)
sample_id <- sub("\\.raw_1\\.fastq\\.gz$", "", basename(r1))
r1_dt <- data.table(sample_id = sample_id, r1_source = normalizePath(r1))
r2_dt <- data.table(
  sample_id = sub("\\.raw_2\\.fastq\\.gz$", "", basename(r2)),
  r2_source = normalizePath(r2)
)
inv <- merge(r1_dt, r2_dt, by = "sample_id", all = TRUE)
if (length(subset_ids) > 0) inv <- inv[sample_id %in% subset_ids]
setorder(inv, sample_id)

if (anyNA(inv$r1_source) || anyNA(inv$r2_source)) {
  stop("Unpaired raw reads detected. Check inventory before continuing.")
}

md5 <- fread(md5_file, header = FALSE, col.names = c("md5", "relative_path"))
md5 <- md5[grepl("\\.raw_[12]\\.fastq\\.gz$", relative_path)]
md5[, sample_id := dirname(relative_path)]
md5[, read := ifelse(grepl("\\.raw_1\\.fastq\\.gz$", relative_path), "R1", "R2")]
md5_wide <- dcast(md5, sample_id ~ read, value.var = "md5")
setnames(md5_wide, old = intersect(names(md5_wide), c("R1", "R2")), new = c("r1_md5_expected", "r2_md5_expected")[seq_along(intersect(names(md5_wide), c("R1", "R2")))])
inv <- merge(inv, md5_wide, by = "sample_id", all.x = TRUE)

info <- rbindlist(lapply(inv$sample_id, function(id) as.data.table(c(sample_id = id, sample_info_from_id(id)))), fill = TRUE)
inv <- merge(inv, info, by = "sample_id", all.x = TRUE)
inv[, `:=`(
  r1_size_bytes = file.info(r1_source)$size,
  r2_size_bytes = file.info(r2_source)$size
)]

inv[, sample_dir := file.path(symlink_root, sample_id)]
inv[, r1_symlink := file.path(sample_dir, basename(r1_source))]
inv[, r2_symlink := file.path(sample_dir, basename(r2_source))]
for (i in seq_len(nrow(inv))) {
  ensure_dir(inv$sample_dir[i])
  if (file.exists(inv$r1_symlink[i]) || nzchar(Sys.readlink(inv$r1_symlink[i]))) unlink(inv$r1_symlink[i])
  if (file.exists(inv$r2_symlink[i]) || nzchar(Sys.readlink(inv$r2_symlink[i]))) unlink(inv$r2_symlink[i])
  file.symlink(inv$r1_source[i], inv$r1_symlink[i])
  file.symlink(inv$r2_source[i], inv$r2_symlink[i])
}

md5_check <- inv[, .(
  sample_id,
  read = c("R1", "R2"),
  file = c(r1_source, r2_source),
  expected_md5 = c(r1_md5_expected, r2_md5_expected)
), by = sample_id]
md5_check[, observed_md5 := unname(tools::md5sum(file))]
md5_check[, status := fifelse(expected_md5 == observed_md5, "pass", "fail")]

sample_metadata <- inv[, .(
  sample_id,
  sample_type,
  include_qc,
  include_biological_analysis,
  diet_code,
  diet_label,
  hydrolysate_percent,
  fishmeal_percent,
  timepoint_day,
  tank,
  fish,
  group,
  forward_raw = r1_source,
  reverse_raw = r2_source,
  forward_symlink = r1_symlink,
  reverse_symlink = r2_symlink
)]

expected <- data.table(sample_id = expected_sample_ids())
if (is_subset_run) expected <- expected[sample_id %in% subset_ids]
expected[, expected := "yes"]
expected_info <- rbindlist(lapply(expected$sample_id, function(id) as.data.table(c(sample_id = id, sample_info_from_id(id)))), fill = TRUE)
expected <- merge(expected, expected_info, by = "sample_id", all.x = TRUE)
expected[, observed_raw_pair := fifelse(sample_id %in% inv$sample_id, "yes", "no")]
missing <- expected[observed_raw_pair == "no"]

dict <- data.table(
  column = names(sample_metadata),
  description = c(
    "Unique sample identifier.",
    "Sample type inferred from sample id.",
    "Whether sample should be retained for QC summaries.",
    "Whether sample is part of biological downstream analysis.",
    "Diet code X/Y/Z where applicable.",
    "Human-readable diet label.",
    "Hydrolysate percentage where applicable.",
    "Fishmeal percentage where applicable.",
    "Sampling day.",
    "Tank number inferred from sample id.",
    "Fish number inferred from sample id.",
    "Diet-timepoint group.",
    "Original raw R1 FASTQ path.",
    "Original raw R2 FASTQ path.",
    "Project-local symlink to raw R1.",
    "Project-local symlink to raw R2."
  )
)

if (!is_subset_run) {
  write_tsv(inv[, !c("sample_dir")], file.path(metadata_dir, "raw_read_inventory.tsv"))
  write_tsv(md5_check, file.path(metadata_dir, "raw_read_md5_check.tsv"))
  write_tsv(sample_metadata, file.path(metadata_dir, "sample_metadata.tsv"))
  write_tsv(expected, file.path(metadata_dir, "expected_samples.tsv"))
  write_tsv(missing, file.path(metadata_dir, "missing_expected_samples.tsv"))
}
write_tsv(inv[, !c("sample_dir")], file.path(metadata_dir, paste0("raw_read_inventory_", args$run_id, ".tsv")))
write_tsv(md5_check, file.path(metadata_dir, paste0("raw_read_md5_check_", args$run_id, ".tsv")))
write_tsv(sample_metadata, file.path(metadata_dir, paste0("sample_metadata_", args$run_id, ".tsv")))
write_tsv(expected, file.path(metadata_dir, paste0("expected_samples_", args$run_id, ".tsv")))
write_tsv(missing, file.path(metadata_dir, paste0("missing_expected_samples_", args$run_id, ".tsv")))
write_tsv(dict, file.path(metadata_dir, "sample_metadata_dictionary.tsv"))

log_message("Raw paired libraries:", nrow(inv))
log_message("MD5 pass:", sum(md5_check$status == "pass"), "fail:", sum(md5_check$status == "fail"))
if (nrow(missing) > 0) log_message("Missing expected samples:", paste(missing$sample_id, collapse = ", "))
if (any(md5_check$status != "pass")) stop("At least one raw FASTQ MD5 check failed.")

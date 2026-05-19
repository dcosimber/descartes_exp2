#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]])), "lib_pipeline.R"))

args <- parse_cli(list(run_id = "run_01_default"))
project <- script_project_root()
run_id <- args$run_id

downstream_sample_info <- function(sample_id) {
  if (grepl("^M[0-9]+$", sample_id)) {
    return(list(
      sample_id_clean = sprintf("Mock_M%02d", as.integer(sub("^M", "", sample_id))),
      diet_code = NA_character_, diet_label = "mock_community",
      treatment = "Mock", treatment_label = "mock_community",
      hydrolysate_percent = NA_integer_, timepoint_label = NA_character_,
      group = "Mock"
    ))
  }

  feed <- grepl("^[XYZ]3$", sample_id)
  m <- regexec("^([XYZ])([0-9]+)\\.([0-9]+)t([0-9]+)$", sample_id)
  parts <- regmatches(sample_id, m)[[1]]
  if (!feed && length(parts) != 5) {
    return(list(
      sample_id_clean = sample_id,
      diet_code = NA_character_, diet_label = NA_character_,
      treatment = NA_character_, treatment_label = NA_character_,
      hydrolysate_percent = NA_integer_, timepoint_label = NA_character_,
      group = "Unknown"
    ))
  }

  code <- if (feed) substr(sample_id, 1, 1) else parts[[2]]
  treatment <- unname(c(X = "Ctrl", Y = "BL15", Z = "BL30")[code])
  label <- unname(c(
    Ctrl = "control",
    BL15 = "blue_whiting_hydrolysate_15",
    BL30 = "blue_whiting_hydrolysate_30"
  )[treatment])
  hydrolysate <- unname(c(Ctrl = 0L, BL15 = 15L, BL30 = 30L)[treatment])

  if (feed) {
    return(list(
      sample_id_clean = paste0(treatment, "_Feed"),
      diet_code = treatment, diet_label = paste0(label, "_feed"),
      treatment = treatment, treatment_label = paste0(label, "_feed"),
      hydrolysate_percent = hydrolysate, timepoint_label = NA_character_,
      group = paste0(treatment, "_Feed")
    ))
  }

  tank <- as.integer(parts[[3]])
  fish <- as.integer(parts[[4]])
  timepoint <- as.integer(parts[[5]])
  timepoint_label <- sprintf("D%02d", timepoint)
  list(
    sample_id_clean = sprintf("%s_%s_Tk%d_F%02d", treatment, timepoint_label, tank, fish),
    diet_code = treatment, diet_label = label,
    treatment = treatment, treatment_label = label,
    hydrolysate_percent = hydrolysate, timepoint_label = timepoint_label,
    group = paste0(treatment, "_", timepoint_label)
  )
}

export_dir <- run_dir(project, "07_exports", run_id)
outdir <- run_dir(project, "09_downstream_ready", run_id)
if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
ensure_dir(outdir)

asv_table <- read_tsv(file.path(export_dir, "asv_table.tsv"))
sample_metadata <- read_tsv(file.path(export_dir, "sample_metadata.tsv"))
sequence_map <- read_tsv(file.path(export_dir, "asv_sequence_map.tsv"))
track_reads <- read_tsv(file.path(export_dir, "track_reads.tsv"))

recoded_info <- rbindlist(lapply(sample_metadata$sample_id, function(id) {
  as.data.table(c(sample_id = id, downstream_sample_info(id)))
}), fill = TRUE)
for (col in setdiff(names(recoded_info), "sample_id")) {
  sample_metadata[[col]] <- recoded_info[[col]][match(sample_metadata$sample_id, recoded_info$sample_id)]
}
if (anyNA(sample_metadata$sample_id_clean) || anyDuplicated(sample_metadata$sample_id_clean)) {
  stop("Clean sample identifiers are missing or duplicated. Check sample metadata before downstream export.")
}
setcolorder(sample_metadata, c("sample_id", "sample_id_clean", setdiff(names(sample_metadata), c("sample_id", "sample_id_clean"))))

write_tsv(asv_table, file.path(outdir, "otu_table_taxa_are_rows.tsv"))
counts <- asv_table[, -"asv_id"]
sample_rows <- as.data.table(t(as.matrix(counts)))
setnames(sample_rows, asv_table$asv_id)
sample_rows[, sample_id := names(counts)]
setcolorder(sample_rows, "sample_id")
write_tsv(sample_rows, file.path(outdir, "otu_table_samples_are_rows.tsv"))
write_tsv(sample_metadata, file.path(outdir, "sample_data.tsv"))
write_tsv(sequence_map, file.path(outdir, "asv_sequence_map.tsv"))
invisible(file.copy(file.path(export_dir, "asv_sequences.fasta"), file.path(outdir, "rep_seqs.fasta"), overwrite = TRUE))
invisible(file.copy(file.path(export_dir, "track_reads.tsv"), file.path(outdir, "preprocessing_track_reads.tsv"), overwrite = TRUE))

tax_files <- list.files(export_dir, pattern = "^taxonomy(_[A-Za-z0-9_.-]+)?\\.tsv$", full.names = TRUE)
tax_files <- tax_files[basename(tax_files) != "taxonomy_database_summary.tsv"]
tax_files <- tax_files[!grepl("^taxonomy_bootstraps", basename(tax_files))]
for (tax_file in tax_files) {
  db_label <- sub("^taxonomy_?", "", tools::file_path_sans_ext(basename(tax_file)))
  if (db_label == "") db_label <- "primary"
  tax <- read_tsv(tax_file)
  tax_cols <- setdiff(names(tax), c("sequence"))
  if ("asv_id" %in% names(tax)) {
    tax_out <- tax[, c("asv_id", setdiff(tax_cols, "asv_id")), with = FALSE]
  } else {
    tax_out <- merge(sequence_map[, .(asv_id, sequence)], tax, by = "sequence", all.x = TRUE)
    tax_out <- tax_out[, c("asv_id", setdiff(names(tax_out), c("asv_id", "sequence"))), with = FALSE]
  }
  write_tsv(tax_out, file.path(outdir, paste0("tax_table_", db_label, ".tsv")))
}

notes <- c(
  paste0("# Downstream-ready inputs for ", run_id),
  "",
  "- `otu_table_taxa_are_rows.tsv`: ASVs in rows, samples in columns; use `taxa_are_rows=TRUE` in phyloseq.",
  "- `otu_table_samples_are_rows.tsv`: samples in rows, ASVs in columns for ordination/statistical tools that prefer sample-wise matrices.",
  "- `sample_data.tsv`: complete sample metadata.",
  "- `tax_table_primary.tsv` or `tax_table_<database>.tsv`: taxonomy tables keyed by `asv_id`.",
  "- `rep_seqs.fasta`: ASV representative sequences with ASV IDs.",
  "- `preprocessing_track_reads.tsv`: read retention/QC by sample.",
  "- `sample_data.tsv` keeps original technical `sample_id` values and adds clean downstream metadata fields such as `sample_id_clean`, `diet_code`, `timepoint_label`, and `group`.",
  "",
  "No downstream normalization, filtering, rarefaction, compositional transform, or phyloseq object construction is applied here."
)
writeLines(notes, file.path(outdir, "README.md"))
log_message("Downstream-ready inputs complete:", outdir)

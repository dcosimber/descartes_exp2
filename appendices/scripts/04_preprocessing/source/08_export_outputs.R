#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(sessioninfo)
})

source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]])), "lib_pipeline.R"))

args <- parse_cli(list(run_id = "run_01_default"))
project <- script_project_root()
run_id <- args$run_id
params <- load_run_params(project, run_id)

metadata <- read_run_metadata(project, run_id)
cutadapt_stats <- read_tsv(file.path(run_dir(project, "01_cutadapt", run_id), "cutadapt_stats.tsv"))
filter_stats <- read_tsv(file.path(run_dir(project, "02_filter_trim", run_id), "dada2_filter_stats.tsv"))
merge_stats <- read_tsv(file.path(run_dir(project, "05_chimera", run_id), "dada2_merge_chimera_stats.tsv"))
seqtab_nochim <- readRDS(file.path(run_dir(project, "05_chimera", run_id), "seqtab_nochim.rds"))
taxonomy <- read_tsv(file.path(run_dir(project, "06_taxonomy", run_id), "taxonomy.tsv"))
boot <- read_tsv(file.path(run_dir(project, "06_taxonomy", run_id), "taxonomy_bootstraps.tsv"))
if (anyDuplicated(names(taxonomy))) names(taxonomy) <- make.unique(names(taxonomy), sep = "_")
if (anyDuplicated(names(boot))) names(boot) <- make.unique(names(boot), sep = "_")
taxonomy_summary_file <- file.path(run_dir(project, "06_taxonomy", run_id), "taxonomy_database_summary.tsv")
taxonomy_summary <- if (file.exists(taxonomy_summary_file)) read_tsv(taxonomy_summary_file) else data.table()

outdir <- run_dir(project, "07_exports", run_id)
report_dir <- file.path(project, "reports", "preprocessing_summary")
ensure_dir(outdir, report_dir, file.path(project, "logs", "versions"))

seqs <- colnames(seqtab_nochim)
asv_ids <- sprintf("ASV%05d", seq_along(seqs))
sequence_map <- data.table(asv_id = asv_ids, sequence = seqs, length = nchar(seqs), total_abundance = colSums(seqtab_nochim))

asv_table <- as.data.table(t(seqtab_nochim))
asv_table[, asv_id := asv_ids]
setcolorder(asv_table, "asv_id")
write_tsv(asv_table, file.path(outdir, "asv_table.tsv"))
write_tsv(sequence_map, file.path(outdir, "asv_sequence_map.tsv"))

fasta_lines <- as.vector(rbind(paste0(">", sequence_map$asv_id), sequence_map$sequence))
writeLines(fasta_lines, file.path(outdir, "asv_sequences.fasta"))

taxonomy <- merge(sequence_map[, .(asv_id, sequence)], taxonomy, by = "sequence", all.x = TRUE)
boot <- merge(sequence_map[, .(asv_id, sequence)], boot, by = "sequence", all.x = TRUE)
write_tsv(taxonomy, file.path(outdir, "taxonomy.tsv"))
write_tsv(boot, file.path(outdir, "taxonomy_bootstraps.tsv"))
if (nrow(taxonomy_summary) > 0) write_tsv(taxonomy_summary, file.path(outdir, "taxonomy_database_summary.tsv"))

tax_run_dir <- run_dir(project, "06_taxonomy", run_id)
tax_db_dirs <- list.dirs(tax_run_dir, recursive = FALSE, full.names = TRUE)
for (db_dir in tax_db_dirs) {
  db_id <- basename(db_dir)
  db_tax_file <- file.path(db_dir, "taxonomy.tsv")
  db_boot_file <- file.path(db_dir, "taxonomy_bootstraps.tsv")
  if (file.exists(db_tax_file)) {
    db_tax <- read_tsv(db_tax_file)
    if (anyDuplicated(names(db_tax))) names(db_tax) <- make.unique(names(db_tax), sep = "_")
    db_tax <- merge(sequence_map[, .(asv_id, sequence)], db_tax, by = "sequence", all.x = TRUE)
    write_tsv(db_tax, file.path(outdir, paste0("taxonomy_", db_id, ".tsv")))
  }
  if (file.exists(db_boot_file)) {
    db_boot <- read_tsv(db_boot_file)
    if (anyDuplicated(names(db_boot))) names(db_boot) <- make.unique(names(db_boot), sep = "_")
    db_boot <- merge(sequence_map[, .(asv_id, sequence)], db_boot, by = "sequence", all.x = TRUE)
    write_tsv(db_boot, file.path(outdir, paste0("taxonomy_bootstraps_", db_id, ".tsv")))
  }
}

track <- merge(metadata[, .(sample_id, sample_type, include_biological_analysis, diet_code, timepoint_day, group)], cutadapt_stats, by = "sample_id", all.x = TRUE)
track <- merge(track, filter_stats[, .(sample_id, reads.in, reads.out)], by = "sample_id", all.x = TRUE)
track <- merge(track, merge_stats, by = "sample_id", all.x = TRUE, suffixes = c("", "_dada2"))
track[, `:=`(
  cutadapt_retention_pct = pairs_written / pairs_processed * 100,
  raw_to_nonchim_pct = nonchim / pairs_processed * 100
)]
write_tsv(track, file.path(outdir, "track_reads.tsv"))

write_tsv(metadata, file.path(outdir, "sample_metadata.tsv"))
write_tsv(cutadapt_stats, file.path(outdir, "cutadapt_stats.tsv"))
write_tsv(filter_stats, file.path(outdir, "dada2_filter_stats.tsv"))
write_tsv(merge_stats, file.path(outdir, "dada2_merge_chimera_stats.tsv"))
write_yaml(params, file.path(outdir, "run_parameters.yaml"))
writeLines(capture.output(sessioninfo::session_info()), file.path(outdir, "session_info.txt"))
writeLines(capture.output(sessioninfo::session_info()), file.path(project, "logs", "versions", paste0(run_id, "_session_info.txt")))

depths <- params$qc$depth_thresholds
bio <- track[include_biological_analysis == "yes"]
threshold_lines <- vapply(depths, function(d) sprintf("- Biological samples below %s reads: %s/%s", d, sum(bio$nonchim < d, na.rm = TRUE), nrow(bio)), character(1))

rank_cols <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), names(taxonomy))
tax_lines <- character()
if (nrow(taxonomy_summary) > 0) {
  tax_lines <- c(tax_lines, capture.output(print(taxonomy_summary)))
} else if (length(rank_cols) > 0) {
  tax_lines <- c(tax_lines, sprintf("- ASVs with %s: %s", rank_cols, colSums(!is.na(taxonomy[, ..rank_cols]) & taxonomy[, ..rank_cols] != "")))
}

report <- c(
  paste0("# ", run_id, " QC report"),
  "",
  "## Summary",
  "",
  sprintf("- Samples: %s", nrow(track)),
  sprintf("- Biological samples: %s", nrow(bio)),
  sprintf("- ASVs: %s", nrow(sequence_map)),
  sprintf("- Non-chimeric reads: %s", sum(track$nonchim, na.rm = TRUE)),
  sprintf("- Median non-chimeric reads/sample: %.0f", median(track$nonchim, na.rm = TRUE)),
  sprintf("- Median cutadapt retention: %.2f%%", median(track$cutadapt_retention_pct, na.rm = TRUE)),
  sprintf("- Median filter retention: %.2f%%", median(track$filter_retention_pct, na.rm = TRUE)),
  sprintf("- Median non-chimeric/input retention: %.2f%%", median(track$nonchim_retention_pct, na.rm = TRUE)),
  "",
  "## Depth thresholds",
  "",
  threshold_lines,
  "",
  "## Sequence lengths",
  "",
  sprintf("- Min: %s", min(sequence_map$length)),
  sprintf("- Median: %.0f", median(sequence_map$length)),
  sprintf("- Max: %s", max(sequence_map$length)),
  "",
  "## Taxonomy",
  "",
  tax_lines,
  "",
  "## Files",
  "",
  "- `asv_table.tsv`",
  "- `asv_sequences.fasta`",
  "- `taxonomy.tsv`",
  "- `taxonomy_<database>.tsv` when multiple databases are configured",
  "- `track_reads.tsv`"
)
writeLines(report, file.path(report_dir, paste0(run_id, "_qc.md")))
log_message("Exports and QC report complete:", outdir)

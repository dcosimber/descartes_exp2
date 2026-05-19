#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dada2)
})

source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]])), "lib_pipeline.R"))

args <- parse_cli(list(run_id = "run_01_default", only_database = ""))
project <- script_project_root()
run_id <- args$run_id
params <- load_run_params(project, run_id)
tax_params <- params$taxonomy

download_if_missing <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 0) return(invisible(TRUE))
  log_message("Downloading", basename(dest))
  download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
}

taxonomy_databases <- function(tax_params) {
  if (!is.null(tax_params$databases)) return(tax_params$databases)
  list(silva_138_2 = list(
    label = tax_params$database %||% "SILVA",
    version = tax_params$version %||% "unknown",
    reference_dir = "silva_138_2_dada2",
    training_set_filename = "silva_nr99_v138.2_toSpecies_trainset.fa.gz",
    training_set_url = tax_params$training_set_url,
    species_set_filename = "silva_v138.2_assignSpecies.fa.gz",
    species_set_url = tax_params$species_set_url,
    min_boot = tax_params$min_boot %||% 50,
    try_rc = tax_params$try_rc %||% FALSE,
    add_species = TRUE
  ))
}

seqtab_nochim <- readRDS(file.path(run_dir(project, "05_chimera", run_id), "seqtab_nochim.rds"))
outdir <- run_dir(project, "06_taxonomy", run_id)
ensure_dir(outdir, file.path(project, "logs", "commands"))
cmd_log <- file.path(project, "logs", "commands", paste0("07_assign_taxonomy_", run_id, ".commands.R"))
writeLines(capture.output(str(tax_params)), cmd_log)

dbs <- taxonomy_databases(tax_params)
primary_db <- tax_params$primary %||% names(dbs)[[1]]
if (!is.null(args$only_database) && !isTRUE(args$only_database) && !is.na(args$only_database) && nzchar(args$only_database)) {
  dbs <- dbs[names(dbs) %in% args$only_database]
  if (length(dbs) == 0) stop("Requested taxonomy database not configured: ", args$only_database)
}
summary_rows <- list()

for (db_id in names(dbs)) {
  db <- dbs[[db_id]]
  ref_dir <- file.path(project, "data", "reference_databases", db$reference_dir %||% db_id)
  db_outdir <- file.path(outdir, db_id)
  ensure_dir(ref_dir, db_outdir)
  train_fa <- file.path(ref_dir, db$training_set_filename)
  download_if_missing(db$training_set_url, train_fa)

  species_fa <- NA_character_
  if (isTRUE(db$add_species) && !is.null(db$species_set_url) && !is.null(db$species_set_filename)) {
    species_fa <- file.path(ref_dir, db$species_set_filename)
    download_if_missing(db$species_set_url, species_fa)
  }

  min_boot <- as.integer(db$min_boot %||% tax_params$min_boot %||% 50)
  try_rc <- isTRUE(db$try_rc %||% tax_params$try_rc)
  multithread <- if (!is.null(db$multithread)) isTRUE(db$multithread) else TRUE
  tax_levels <- db$tax_levels

  log_message("Assigning taxonomy with", db_id, basename(train_fa))
  assign_args <- list(
    seqs = seqtab_nochim,
    refFasta = train_fa,
    minBoot = min_boot,
    tryRC = try_rc,
    multithread = multithread,
    outputBootstraps = TRUE
  )
  if (!is.null(tax_levels)) assign_args$taxLevels <- unlist(tax_levels)
  taxa <- do.call(assignTaxonomy, assign_args)
  tax_table <- taxa$tax
  boot_table <- taxa$boot

  tax_species <- NULL
  if (!is.na(species_fa) && file.exists(species_fa) && file.info(species_fa)$size > 0) {
    log_message("Adding species with", db_id, basename(species_fa))
    tax_species <- addSpecies(tax_table, species_fa, tryRC = try_rc, allowMultiple = TRUE)
  }

  saveRDS(tax_table, file.path(db_outdir, paste0(db_id, "_assignTaxonomy.rds")))
  saveRDS(boot_table, file.path(db_outdir, paste0(db_id, "_bootstraps.rds")))
  if (!is.null(tax_species)) saveRDS(tax_species, file.path(db_outdir, paste0(db_id, "_addSpecies.rds")))

  tax_final <- if (!is.null(tax_species)) tax_species else tax_table
  tax_export <- as.data.table(tax_final, keep.rownames = "sequence")
  if (anyDuplicated(names(tax_export))) {
    names(tax_export) <- make.unique(names(tax_export), sep = "_")
    names(tax_export) <- sub("^Species_1$", "Species_addSpecies", names(tax_export))
  }
  boot_export <- as.data.table(boot_table, keep.rownames = "sequence")
  write_tsv(tax_export, file.path(db_outdir, "taxonomy.tsv"))
  write_tsv(boot_export, file.path(db_outdir, "taxonomy_bootstraps.tsv"))

  if (identical(db_id, primary_db)) {
    write_tsv(tax_export, file.path(outdir, "taxonomy.tsv"))
    write_tsv(boot_export, file.path(outdir, "taxonomy_bootstraps.tsv"))
    saveRDS(tax_table, file.path(outdir, "taxonomy_assignTaxonomy.rds"))
    saveRDS(boot_table, file.path(outdir, "taxonomy_bootstraps.rds"))
    if (!is.null(tax_species)) saveRDS(tax_species, file.path(outdir, "taxonomy_addSpecies.rds"))
  }

  rank_cols <- intersect(c("Kingdom", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species_addSpecies"), names(tax_export))
  assigned <- if (length(rank_cols) > 0) {
    tax_export[, lapply(.SD, function(x) sum(!is.na(x) & x != "")), .SDcols = rank_cols]
  } else data.table()
  summary_rows[[db_id]] <- data.table(
    database_id = db_id,
    label = db$label %||% db_id,
    version = db$version %||% NA_character_,
    min_boot = min_boot,
    try_rc = try_rc,
    multithread = multithread,
    asvs = nrow(tax_export),
    assigned
  )
  log_message("Taxonomy complete for", db_id, nrow(tax_export), "ASVs.")
}

new_summary <- rbindlist(summary_rows, fill = TRUE)
old_summary_file <- file.path(outdir, "taxonomy_database_summary.tsv")
if (file.exists(old_summary_file)) {
  old_summary <- read_tsv(old_summary_file)
  old_summary <- old_summary[!database_id %in% new_summary$database_id]
  new_summary <- rbindlist(list(old_summary, new_summary), fill = TRUE)
}
write_tsv(new_summary, old_summary_file)
log_message("All taxonomy databases complete for", run_id)

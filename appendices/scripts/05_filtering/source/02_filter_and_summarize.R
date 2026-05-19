#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(phyloseq)
})

source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]])), "lib_downstream.R"))

args <- parse_cli(list(run_id = NULL, taxonomy_database = NULL))
project <- script_project_root()
params <- load_params(project)
run_id <- args$run_id %||% params$project$final_preprocessing_run
tax_db <- args$taxonomy_database %||% params$project$taxonomy_database
objects_dir <- stage_objects_dir(project, "02_filtering")
tables_dir <- stage_tables_dir(project, "02_filtering")
reports_dir <- stage_dir(project, "02_filtering", "reports")
ensure_dir(objects_dir, tables_dir, reports_dir)

ps <- read_object(project, "01_import", "ps_all", tax_db, run_id)
steps <- list(raw_all_samples = ps)

tax <- as.data.frame(tax_table(ps))
tax_dt <- as.data.table(tax, keep.rownames = "taxon_id")
seq_map_file <- file.path(downstream_ready_dir(project, params, run_id), "asv_sequence_map.tsv")
seq_map <- if (file.exists(seq_map_file)) read_tsv(seq_map_file) else data.table(asv_id = taxa_names(ps))
if ("asv_id" %in% names(seq_map)) setnames(seq_map, "asv_id", "taxon_id")

keep <- rep(TRUE, ntaxa(ps))
names(keep) <- taxa_names(ps)
removal_flags <- data.table(taxon_id = taxa_names(ps))
removal_flags[, `:=`(
  remove_eukaryota_unassigned = FALSE,
  remove_mitochondria_chloroplast = FALSE,
  remove_unclassified_all_ranks = FALSE,
  remove_kingdom_only = FALSE
)]
if (isTRUE(params$filtering$remove_eukaryota) && "Kingdom" %in% names(tax)) {
  flag <- tax$Kingdom %in% c("Eukaryota", "Unassigned")
  flag[is.na(flag)] <- FALSE
  removal_flags[, remove_eukaryota_unassigned := flag]
  keep <- keep & !flag
}
if (isTRUE(params$filtering$remove_mitochondria_chloroplast)) {
  bad <- rep(FALSE, ntaxa(ps))
  for (col in intersect(c("Order", "Family", "Genus"), names(tax))) {
    bad <- bad | grepl("chloroplast|mitochond", tax[[col]], ignore.case = TRUE)
  }
  bad[is.na(bad)] <- FALSE
  removal_flags[, remove_mitochondria_chloroplast := bad]
  keep <- keep & !bad
}
if (isTRUE(params$filtering$remove_unclassified_all_ranks)) {
  ranks_all <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), names(tax))
  if (length(ranks_all) > 0) {
    no_assignment <- rowSums(!is.na(tax[, ranks_all, drop = FALSE]) & tax[, ranks_all, drop = FALSE] != "") == 0
    no_assignment[is.na(no_assignment)] <- FALSE
    removal_flags[, remove_unclassified_all_ranks := no_assignment]
    keep <- keep & !no_assignment
  }
}
if (isTRUE(params$filtering$remove_kingdom_only)) {
  ranks_below <- intersect(c("Phylum", "Class", "Order", "Family", "Genus", "Species"), names(tax))
  if (length(ranks_below) > 0 && "Kingdom" %in% names(tax)) {
    kingdom_only <- !is.na(tax$Kingdom) & rowSums(!is.na(tax[, ranks_below, drop = FALSE]) & tax[, ranks_below, drop = FALSE] != "") == 0
    kingdom_only[is.na(kingdom_only)] <- FALSE
    removal_flags[, remove_kingdom_only := kingdom_only]
    keep <- keep & !kingdom_only
  }
}
ps_tax_filtered <- prune_taxa(taxa_names(ps)[keep], ps)
ps_tax_filtered <- prune_taxa(taxa_sums(ps_tax_filtered) > 0, ps_tax_filtered)
steps$taxonomic_filter <- ps_tax_filtered

otu <- as(otu_table(ps_tax_filtered), "matrix")
if (!taxa_are_rows(ps_tax_filtered)) otu <- t(otu)
sample_groups <- sample_data_dt(ps_tax_filtered)
sample_groups[, prevalence_group := fifelse(
  sample_type == params$design$biological_sample_type,
  as.character(diet_time),
  fifelse(
    sample_type == "Feed",
    paste0("Feed_", as.character(diet_code)),
    fifelse(sample_type == "Mock", "Mock", as.character(sample_type))
  )
)]
sample_groups[is.na(prevalence_group) | prevalence_group == "", prevalence_group := sample_type]
sample_groups[, prevalence_group := factor(prevalence_group)]
prevalence_groups <- sample_groups[, .(n_samples = .N), by = prevalence_group][order(prevalence_group)]
prevalence_groups[, `:=`(
  min_prevalence_fraction = params$filtering$min_prevalence_fraction_within_diet_time,
  min_prevalence_absolute_config = params$filtering$min_prevalence_absolute,
  min_prevalence_required = pmin(
    n_samples,
    pmax(
      params$filtering$min_prevalence_absolute,
      ceiling(params$filtering$min_prevalence_fraction_within_diet_time * n_samples)
    )
  )
)]
prev_keep <- rep(FALSE, nrow(otu))
names(prev_keep) <- rownames(otu)
for (grp in na.omit(unique(sample_groups$prevalence_group))) {
  ids <- sample_groups[prevalence_group == grp, sample_id]
  min_prev <- prevalence_groups[prevalence_group == grp, min_prevalence_required]
  prev_keep <- prev_keep | rowSums(otu[, ids, drop = FALSE] > 0) >= min_prev
}
ps_prev <- prune_taxa(names(prev_keep)[prev_keep], ps_tax_filtered)
ps_prev <- prune_taxa(taxa_sums(ps_prev) > 0, ps_prev)
steps$prevalence_filter <- ps_prev

rel_mean <- taxa_sums(transform_sample_counts(ps_prev, function(x) x / sum(x))) / nsamples(ps_prev)
ps_final <- prune_taxa(names(rel_mean)[rel_mean >= params$filtering$min_mean_relative_abundance], ps_prev)
ps_final <- prune_taxa(taxa_sums(ps_final) > 0, ps_final)
steps$mean_abundance_filter <- ps_final

ps_genus <- tax_glom(ps_final, taxrank = "Genus", NArm = FALSE)
ps_family <- tax_glom(ps_final, taxrank = "Family", NArm = FALSE)
ps_phylum <- tax_glom(ps_final, taxrank = "Phylum", NArm = FALSE)

biological_keep <- sample_data(ps_final)$sample_type == params$design$biological_sample_type &
  sample_data(ps_final)$include_biological_analysis == "yes"
ps_biological_final <- prune_samples(biological_keep, ps_final)
ps_biological_final <- prune_taxa(taxa_sums(ps_biological_final) > 0, ps_biological_final)
ps_biological_genus <- tax_glom(ps_biological_final, taxrank = "Genus", NArm = FALSE)
ps_biological_family <- tax_glom(ps_biological_final, taxrank = "Family", NArm = FALSE)
ps_biological_phylum <- tax_glom(ps_biological_final, taxrank = "Phylum", NArm = FALSE)

summarize_ps <- function(name, obj, previous = NULL) {
  reads <- sum(sample_sums(obj))
  taxa <- ntaxa(obj)
  data.table(
    step = name,
    samples = nsamples(obj),
    taxa = taxa,
    reads = reads,
    min_reads = min(sample_sums(obj)),
    median_reads = median(sample_sums(obj)),
    max_reads = max(sample_sums(obj)),
    sparsity = sum(as(otu_table(obj), "matrix") == 0) / length(as(otu_table(obj), "matrix")),
    reads_removed_from_previous = if (is.null(previous)) NA_real_ else sum(sample_sums(previous)) - reads,
    taxa_removed_from_previous = if (is.null(previous)) NA_real_ else ntaxa(previous) - taxa
  )
}
filter_stats <- rbindlist(Map(summarize_ps, names(steps), steps, c(list(NULL), steps[-length(steps)])))
write_tsv(filter_stats, file.path(tables_dir, "filtering_summary.tsv"))

summarize_by_sample_type <- function(step_name, obj) {
  meta <- sample_data_dt(obj)
  dt <- data.table(sample_id = sample_names(obj), reads = sample_sums(obj))
  dt <- merge(dt, meta[, .(sample_id, sample_type)], by = "sample_id")
  out <- dt[, .(
    samples = .N,
    reads = sum(reads),
    min_reads = min(reads),
    median_reads = median(reads),
    max_reads = max(reads)
  ), by = sample_type]
  out[, step := step_name]
  setcolorder(out, c("step", "sample_type"))
  out
}
filter_stats_by_sample_type <- rbindlist(Map(summarize_by_sample_type, names(steps), steps), fill = TRUE)
write_tsv(filter_stats_by_sample_type, file.path(tables_dir, "filtering_summary_by_sample_type.tsv"))

sample_depths <- data.table(sample_id = sample_names(ps_final), final_reads = sample_sums(ps_final))
sample_depths <- merge(sample_data_dt(ps_final), sample_depths, by = "sample_id")
write_tsv(sample_depths, file.path(tables_dir, "final_sample_depths.tsv"))

biological_sample_depths <- data.table(sample_id = sample_names(ps_biological_final), final_reads = sample_sums(ps_biological_final))
biological_sample_depths <- merge(sample_data_dt(ps_biological_final), biological_sample_depths, by = "sample_id")
write_tsv(biological_sample_depths, file.path(tables_dir, "final_biological_sample_depths.tsv"))

sample_group_table <- sample_groups[, .(sample_id, sample_type, diet_code, timepoint, diet_time, prevalence_group)]
write_tsv(sample_group_table, file.path(tables_dir, "sample_filtering_groups.tsv"))
write_tsv(prevalence_groups, file.path(tables_dir, "prevalence_filter_groups.tsv"))

taxa_status <- data.table(taxon_id = taxa_names(ps))
taxa_status[, `:=`(
  raw_reads = as.numeric(taxa_sums(ps)[taxon_id]),
  after_taxonomic_filter = taxon_id %in% taxa_names(ps_tax_filtered),
  after_prevalence_filter = taxon_id %in% taxa_names(ps_prev),
  after_mean_abundance_filter = taxon_id %in% taxa_names(ps_final)
)]
taxa_status <- merge(taxa_status, removal_flags, by = "taxon_id", all.x = TRUE)
taxa_status <- merge(taxa_status, tax_dt, by = "taxon_id", all.x = TRUE)
taxa_status <- merge(taxa_status, seq_map, by = "taxon_id", all.x = TRUE)
taxa_status[, removal_step := fifelse(
  after_mean_abundance_filter,
  "retained_final",
  fifelse(
    !after_taxonomic_filter,
    "taxonomic_filter",
    fifelse(!after_prevalence_filter, "prevalence_filter", "mean_abundance_filter")
  )
)]
taxa_status[, removal_reason := removal_step]
taxa_status[removal_step == "retained_final", removal_reason := "retained"]
taxa_status[removal_step == "taxonomic_filter" & remove_eukaryota_unassigned == TRUE, removal_reason := "eukaryota_or_unassigned"]
taxa_status[removal_step == "taxonomic_filter" & remove_mitochondria_chloroplast == TRUE, removal_reason := "mitochondria_or_chloroplast"]
taxa_status[removal_step == "taxonomic_filter" & remove_unclassified_all_ranks == TRUE, removal_reason := "unclassified_all_ranks"]
taxa_status[removal_step == "taxonomic_filter" & remove_kingdom_only == TRUE, removal_reason := "kingdom_only"]
taxa_status[removal_step == "taxonomic_filter" & remove_eukaryota_unassigned == TRUE & remove_mitochondria_chloroplast == TRUE, removal_reason := "eukaryota_or_unassigned;mitochondria_or_chloroplast"]
write_tsv(taxa_status[order(removal_step, -raw_reads)], file.path(tables_dir, "asv_filtering_status.tsv"))
write_tsv(taxa_status[removal_step != "retained_final"][order(removal_step, -raw_reads)], file.path(tables_dir, "removed_asvs.tsv"))
write_tsv(taxa_status[removal_step == "retained_final"][order(-raw_reads)], file.path(tables_dir, "retained_asvs.tsv"))

asv_counts_by_removal_step <- taxa_status[, .(
  asvs = .N,
  raw_reads = sum(raw_reads)
), by = .(removal_step, removal_reason)][order(removal_step, -raw_reads)]
write_tsv(asv_counts_by_removal_step, file.path(tables_dir, "asv_counts_by_filter_step.tsv"))

tax_ranks_available <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), names(taxa_status))
taxonomic_removal_summary <- taxa_status[, .(
  asvs = .N,
  raw_reads = sum(raw_reads)
), by = c("removal_step", tax_ranks_available)][order(removal_step, -raw_reads)]
write_tsv(taxonomic_removal_summary, file.path(tables_dir, "taxonomic_summary_by_filter_step.tsv"))

saveRDS(ps_final, object_file(project, "02_filtering", "ps_final", tax_db))
saveRDS(ps_genus, object_file(project, "02_filtering", "ps_genus", tax_db))
saveRDS(ps_family, object_file(project, "02_filtering", "ps_family", tax_db))
saveRDS(ps_phylum, object_file(project, "02_filtering", "ps_phylum", tax_db))
saveRDS(ps_biological_final, object_file(project, "02_filtering", "ps_biological_final", tax_db))
saveRDS(ps_biological_genus, object_file(project, "02_filtering", "ps_biological_genus", tax_db))
saveRDS(ps_biological_family, object_file(project, "02_filtering", "ps_biological_family", tax_db))
saveRDS(ps_biological_phylum, object_file(project, "02_filtering", "ps_biological_phylum", tax_db))
saveRDS(steps, object_file(project, "02_filtering", "filtering_steps", tax_db))

write_xlsx(
  list(
    filtering_summary = filter_stats,
    summary_by_sample_type = filter_stats_by_sample_type,
    prevalence_groups = prevalence_groups,
    sample_groups = sample_group_table,
    sample_depths = sample_depths,
    biological_sample_depths = biological_sample_depths,
    asv_counts_by_step = asv_counts_by_removal_step,
    removed_asvs = taxa_status[removal_step != "retained_final"][order(removal_step, -raw_reads)],
    retained_asvs = taxa_status[removal_step == "retained_final"][order(-raw_reads)]
  ),
  file.path(tables_dir, "filtering_summary.xlsx")
)

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_pct <- function(x) sprintf("%.2f%%", 100 * x)
final_row <- filter_stats[step == "mean_abundance_filter"]
raw_row <- filter_stats[step == "raw_all_samples"]
bio_final_row <- data.table(
  samples = nsamples(ps_biological_final),
  taxa = ntaxa(ps_biological_final),
  reads = sum(sample_sums(ps_biological_final))
)
report_lines <- c(
  "# Filtering report",
  "",
  paste0("Run: `", run_id, "`  "),
  paste0("Base de datos taxonomica: `", tax_db, "`"),
  "",
  "## 1. Objetivo",
  "",
  "Este bloque genera los objetos phyloseq filtrados que alimentan el analisis downstream. El filtrado se aplica ahora a todos los tipos de muestra importados: intestino, mock community y pienso/feed. Ademas, se exporta un subconjunto biologico intestinal para los analisis que deben excluir controles mock y muestras de pienso.",
  "",
  "## 2. Inputs",
  "",
  "- Objeto de entrada: `results/01_import/objects/ps_all_silva_138_2.rds`",
  "- Tablas fuente: tabla OTU, tabla taxonomica SILVA, metadatos estandarizados y trazabilidad de lecturas del export downstream-ready de DADA2.",
  "",
  "## 3. Diseno del filtrado",
  "",
  "- El filtro taxonomico elimina ASVs asignadas a Eukaryota/Unassigned y ASVs clasificadas solo a nivel de reino, segun `config/downstream_params.yaml`. Cloroplasto y mitocondria se conservan para evaluar la composicion del pienso/feed.",
  "- Las ASVs sin ninguna asignacion taxonomica en los rangos disponibles (`Kingdom` a `Species`) se eliminan como `unclassified_all_ranks`; esto no afecta a cloroplasto/mitocondria porque conservan asignacion en rangos superiores.",
  "- El filtro de prevalencia se evalua dentro de grupos definidos por tipo de muestra. Las muestras intestinales se agrupan por `diet_time`; los piensos como `Feed_Ctrl`, `Feed_BL15` y `Feed_BL30`; y los mocks como `Mock`.",
  paste0("- La regla de prevalencia configurada es max(`", params$filtering$min_prevalence_absolute, "` muestras, `", params$filtering$min_prevalence_fraction_within_diet_time, "` del tamano del grupo), limitada al tamano real del grupo para que los grupos de pienso con una sola muestra puedan retener ASVs presentes."),
  paste0("- El filtro de abundancia media retiene ASVs con abundancia relativa media >= `", params$filtering$min_mean_relative_abundance, "`."),
  "",
  "## 4. Resultados principales",
  "",
  paste0("- Objeto bruto con todas las muestras: ", fmt_int(raw_row$samples), " muestras, ", fmt_int(raw_row$taxa), " ASVs y ", fmt_int(raw_row$reads), " lecturas."),
  paste0("- Objeto final con todas las muestras: ", fmt_int(final_row$samples), " muestras, ", fmt_int(final_row$taxa), " ASVs y ", fmt_int(final_row$reads), " lecturas."),
  paste0("- Subconjunto final biologico intestinal: ", fmt_int(bio_final_row$samples), " muestras, ", fmt_int(bio_final_row$taxa), " ASVs y ", fmt_int(bio_final_row$reads), " lecturas."),
  "",
  "### Tabla 1. Resumen del filtrado",
  "",
  "`tables/filtering_summary.tsv` y la pestana `filtering_summary` de `tables/filtering_summary.xlsx` muestran la retencion de muestras, ASVs y lecturas en cada paso.",
  "",
  "### Tabla 2. Grupos de prevalencia",
  "",
  "`tables/prevalence_filter_groups.tsv` documenta los grupos de muestras y el umbral minimo de prevalencia usado en cada grupo. Esto es importante porque cada grupo de pienso contiene una sola muestra, mientras que el mock tiene dos muestras de control.",
  "",
  "### Tabla 3. ASVs eliminadas y retenidas",
  "",
  "`tables/asv_filtering_status.tsv` contiene una fila por ASV con taxonomia, secuencia cuando esta disponible, recuento de lecturas, flags de retencion y primer paso en el que se elimina. `tables/removed_asvs.tsv` y `tables/retained_asvs.tsv` separan esta tabla de auditoria para facilitar la revision.",
  "",
  "### Tabla 4. Resumen taxonomico por paso de filtrado",
  "",
  "`tables/taxonomic_summary_by_filter_step.tsv` resume que grupos taxonomicos se eliminan o retienen en cada paso. Es la primera tabla que conviene revisar si un taxon relevante desaparece tras el filtrado.",
  "",
  "## 5. Objetos exportados",
  "",
  "- Objetos filtrados globales: `ps_final`, `ps_genus`, `ps_family`, `ps_phylum`.",
  "- Objetos filtrados biologicos intestinales: `ps_biological_final`, `ps_biological_genus`, `ps_biological_family`, `ps_biological_phylum`.",
  "- Trayectoria completa del filtrado: `filtering_steps`.",
  "",
  "## 6. Notas y limitaciones",
  "",
  "- Las muestras mock y feed estan incluidas en la auditoria global del filtrado, pero los analisis biologicos de diversidad alfa, beta y abundancia diferencial deben usar los objetos `ps_biological_*`, salvo que la comparacion incluya explicitamente pienso.",
  "- Las etiquetas a nivel especie siguen siendo exploratorias en datos 16S; el filtrado es a nivel ASV y no depende de la confianza a especie.",
  "- Cloroplasto y mitocondria se conservan por defecto en este proyecto porque aportan informacion relevante sobre la composicion del pienso/feed. Si se necesita una matriz estrictamente bacteriana, debe generarse como analisis de sensibilidad separado.",
  ""
)
# The curated report is maintained at results/02_filtering/reports/filtering_report.md.
# Do not emit an additional auto-report, to avoid conflicting or obsolete reports.

log_message("Filtering complete:", run_id, tax_db)

#!/usr/bin/env Rscript

project_root <- normalizePath(file.path(getwd()), mustWork = TRUE)
source_root <- Sys.getenv(
  "DESCARTES_DOWNSTREAM",
  "/mnt/lustre/scratch/nlsas/home/otras/pia/dci/descartes_exp2/custom_downstream_analysis"
)
source_root <- normalizePath(source_root, mustWork = TRUE)

reports <- data.frame(
  chapter = c(
    "01_filtering",
    "02_dataset_stats",
    "03_rarefaction_qc",
    "04_taxonomic_composition",
    "05_alpha_diversity",
    "06_beta_diversity",
    "07_differential_abundance",
    "08_functional_prediction"
  ),
  stage = c(
    "02_filtering",
    "03_dataset_stats",
    "04_rarefaction_qc",
    "05_taxonomic_composition",
    "06_alpha_diversity",
    "07_beta_diversity",
    "08_differential_abundance",
    "09_functional_prediction"
  ),
  report = c(
    "filtering_report.md",
    "dataset_stats_report.md",
    "rarefaction_qc_report.md",
    "taxonomic_composition_report.md",
    "alpha_diversity_report.md",
    "beta_diversity_report.md",
    "differential_abundance_report.md",
    "functional_prediction_report.md"
  ),
  stringsAsFactors = FALSE
)

rewrite_links <- function(text, stage) {
  base <- paste0("../../assets/results/", stage, "/")
  text <- gsub("\\.\\./figures/", paste0(base, "figures/"), text, fixed = FALSE)
  text <- gsub("\\.\\./tables/", paste0(base, "tables/"), text, fixed = FALSE)
  text <- gsub("\\.\\./objects/", paste0(base, "objects/"), text, fixed = FALSE)
  text <- gsub("\\.\\./logs/", paste0(base, "logs/"), text, fixed = FALSE)
  text <- gsub("\\.\\./picrust2/", paste0(base, "picrust2/"), text, fixed = FALSE)
  text <- gsub("../../03_dataset_stats/reports/dataset_stats_report.md", "02_dataset_stats.qmd", text, fixed = TRUE)
  text <- gsub("../03_dataset_stats/reports/dataset_stats_report.md", "02_dataset_stats.qmd", text, fixed = TRUE)
  text <- gsub("\\.svg", ".png", text)
  text <- gsub("\\.pdf", ".png", text)
  text <- gsub("\\.tsv", ".csv", text)
  text <- text[!grepl("\\.xlsx", text, ignore.case = TRUE)]
  text
}

relative_path <- function(files, root) {
  substring(normalizePath(files, mustWork = TRUE), nchar(normalizePath(root, mustWork = TRUE)) + 2L)
}

copy_figures_as_png <- function(from, to) {
  if (!dir.exists(from)) return(data.frame())
  if (dir.exists(to)) unlink(to, recursive = TRUE, force = TRUE)
  dir.create(to, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(from, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir == FALSE]
  ext <- tolower(tools::file_ext(files))
  files <- files[ext %in% c("svg", "png", "jpg", "jpeg")]
  if (length(files) == 0) return(data.frame())

  rel <- relative_path(files, from)
  stems <- sub("\\.[^.]+$", "", rel)
  ext <- tolower(tools::file_ext(files))

  # Prefer vector SVG sources when both SVG and raster versions share a stem.
  keep <- !duplicated(stems) | ext == "svg"
  files <- files[keep]
  rel <- rel[keep]
  stems <- stems[keep]
  ext <- ext[keep]

  manifest <- data.frame(source = files, output = character(length(files)), stringsAsFactors = FALSE)
  for (i in seq_along(files)) {
    out <- file.path(to, paste0(stems[i], ".png"))
    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    if (ext[i] == "svg") {
      status <- system2("rsvg-convert", c("-o", out, files[i]), stdout = TRUE, stderr = TRUE)
      if (!file.exists(out)) {
        stop("Could not convert SVG to PNG: ", files[i], "\n", paste(status, collapse = "\n"))
      }
    } else {
      file.copy(files[i], out, overwrite = TRUE, copy.date = TRUE)
    }
    manifest$output[i] <- out
  }
  manifest
}

copy_tables_as_csv <- function(from, to) {
  if (!dir.exists(from)) return(data.frame())
  if (dir.exists(to)) unlink(to, recursive = TRUE, force = TRUE)
  dir.create(to, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(from, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir == FALSE]
  ext <- tolower(tools::file_ext(files))
  files <- files[ext %in% c("tsv", "csv")]
  if (length(files) == 0) return(data.frame())

  rel <- relative_path(files, from)
  out_rel <- sub("\\.[^.]+$", ".csv", rel)
  manifest <- data.frame(source = files, output = file.path(to, out_rel), stringsAsFactors = FALSE)

  for (i in seq_along(files)) {
    dir.create(dirname(manifest$output[i]), recursive = TRUE, showWarnings = FALSE)
    sep <- if (tolower(tools::file_ext(files[i])) == "tsv") "\t" else ","
    tab <- utils::read.delim(files[i], sep = sep, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "\"")
    utils::write.csv(tab, manifest$output[i], row.names = FALSE, na = "")
  }
  manifest
}

assets_root <- file.path(project_root, "assets", "results")
generated_root <- file.path(project_root, "chapters", "generated")
if (dir.exists(assets_root)) unlink(assets_root, recursive = TRUE, force = TRUE)
dir.create(assets_root, recursive = TRUE, showWarnings = FALSE)
dir.create(generated_root, recursive = TRUE, showWarnings = FALSE)

manifest <- data.frame(
  chapter = character(),
  stage = character(),
  source_report = character(),
  output_chapter = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(reports))) {
  stage <- reports$stage[i]
  stage_source <- file.path(source_root, "results", stage)
  source_report <- file.path(stage_source, "reports", reports$report[i])
  output_chapter <- file.path(generated_root, paste0(reports$chapter[i], ".qmd"))

  if (!file.exists(source_report)) {
    warning("Missing source report: ", source_report)
    next
  }

  stage_assets <- file.path(assets_root, stage)
  dir.create(stage_assets, recursive = TRUE, showWarnings = FALSE)
  copy_figures_as_png(file.path(stage_source, "figures"), file.path(stage_assets, "figures"))
  copy_tables_as_csv(file.path(stage_source, "tables"), file.path(stage_assets, "tables"))
  copy_figures_as_png(file.path(stage_source, "picrust2", "figures"), file.path(stage_assets, "picrust2", "figures"))
  copy_tables_as_csv(file.path(stage_source, "picrust2", "tables"), file.path(stage_assets, "picrust2", "tables"))

  text <- readLines(source_report, warn = FALSE)
  text <- rewrite_links(text, stage)
  writeLines(text, output_chapter, useBytes = TRUE)

  manifest <- rbind(
    manifest,
    data.frame(
      chapter = reports$chapter[i],
      stage = stage,
      source_report = source_report,
      output_chapter = output_chapter,
      stringsAsFactors = FALSE
    )
  )
}

utils::write.csv(
  manifest,
  file.path(project_root, "assets", "sync_manifest.csv"),
  row.names = FALSE,
  na = ""
)

message("Synced ", nrow(manifest), " report chapters from ", source_root)

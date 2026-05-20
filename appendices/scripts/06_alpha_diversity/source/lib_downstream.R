suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

script_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[startsWith(args, file_arg)][1])
  if (is.na(script_path)) normalizePath(getwd()) else normalizePath(file.path(dirname(script_path), ".."))
}

parse_cli <- function(defaults = list()) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- gsub("-", "_", sub("^--", "", key))
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1
    } else {
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

ensure_dir <- function(...) {
  dirs <- unlist(list(...), use.names = FALSE)
  for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

read_tsv <- function(file) data.table::fread(file, sep = "\t", na.strings = c("NA", ""))
write_tsv <- function(x, file) {
  ensure_dir(dirname(file))
  data.table::fwrite(as.data.table(x), file, sep = "\t", na = "NA", quote = FALSE)
}

write_xlsx <- function(sheets, file) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    warning("Skipping Excel export because openxlsx is not installed: ", file)
    return(invisible(FALSE))
  }
  ensure_dir(dirname(file))
  wb <- openxlsx::createWorkbook()
  for (sheet_name in names(sheets)) {
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, sheets[[sheet_name]])
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  invisible(TRUE)
}

load_params <- function(project = script_project_root()) {
  yaml::read_yaml(file.path(project, "config", "downstream_params.yaml"))
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_message <- function(...) cat(sprintf("[%s] %s\n", timestamp(), paste(..., collapse = " ")))

preprocessing_export_dir <- function(project, params, run_id = NULL) {
  run_id <- run_id %||% params$project$final_preprocessing_run
  normalizePath(file.path(project, params$project$preprocessing_project, "results", "07_exports", run_id), mustWork = TRUE)
}

downstream_ready_dir <- function(project, params, run_id = NULL) {
  run_id <- run_id %||% params$project$final_preprocessing_run
  normalizePath(file.path(project, params$project$preprocessing_project, "results", "09_downstream_ready", run_id), mustWork = TRUE)
}

result_dir <- function(project, ...) {
  file.path(project, "results", ...)
}

report_dir <- function(project, ...) {
  file.path(project, "reports", ...)
}

stage_dir <- function(project, stage, ...) {
  file.path(project, "results", stage, ...)
}

stage_tables_dir <- function(project, stage) stage_dir(project, stage, "tables")
stage_figures_dir <- function(project, stage) stage_dir(project, stage, "figures")
stage_objects_dir <- function(project, stage) stage_dir(project, stage, "objects")
stage_logs_dir <- function(project, stage) stage_dir(project, stage, "logs")

object_file <- function(project, stage, name, tax_db) {
  file.path(stage_objects_dir(project, stage), paste0(name, "_", tax_db, ".rds"))
}

read_object <- function(project, stage, name, tax_db, run_id = NULL) {
  candidates <- c(object_file(project, stage, name, tax_db))
  if (!is.null(run_id) && !is.na(run_id) && nzchar(run_id)) {
    candidates <- c(candidates, file.path(result_dir(project, "objects"), paste0(run_id, "_", tax_db, "_", name, ".rds")))
  }
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("Missing object ", name, " for ", tax_db, ". Checked: ", paste(candidates, collapse = ", "))
  }
  readRDS(existing[[1]])
}

clean_tax_string <- function(x) {
  x <- as.character(x)
  x <- gsub("^[dkpcofgs]__", "", x)
  x <- gsub("^D_[0-9]__", "", x)
  x[x == "" | is.na(x)] <- NA_character_
  x
}

standardize_metadata <- function(meta, params) {
  meta <- as.data.table(meta)
  meta[, diet_code := factor(diet_code, levels = params$design$diet_levels)]
  meta[, timepoint_day := as.integer(timepoint_day)]
  meta[, timepoint_label := fifelse(is.na(timepoint_day), NA_character_, sprintf("D%02d", timepoint_day))]
  meta[, timepoint := factor(timepoint_label, levels = sprintf("D%02d", params$design$timepoint_levels))]
  hydro_levels <- params$design$hydro_group_levels %||% c("Ctrl", "Hydrolysate")
  meta[, hydro_group := fifelse(
    as.character(diet_code) == "Ctrl",
    "Ctrl",
    fifelse(as.character(diet_code) %in% c("BL15", "BL30"), "Hydrolysate", NA_character_)
  )]
  meta[, hydro_group := factor(hydro_group, levels = hydro_levels)]
  diet_labels <- unlist(params$design$diet_labels)
  meta[, diet_label_short := diet_labels[as.character(diet_code)]]
  meta[, diet_time := factor(paste0(as.character(diet_code), "_", timepoint_label),
                             levels = names(unlist(params$colors$diet_time)))]
  meta[, tank := factor(tank)]
  meta[, fish := factor(fish)]
  meta
}

theme_publication <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_line(linewidth = 0.18, color = "grey88"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey65", linewidth = 0.25),
      axis.text = element_text(color = "grey15"),
      legend.key = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey35")
    )
}

diet_colors <- function(params) unlist(params$colors$diet)
hydro_colors <- function(params) {
  cols <- unlist(params$colors$hydro_group)
  if (length(cols) == 0) cols <- c(Ctrl = "#2B6CB0", Hydrolysate = "#C53030")
  cols
}
time_colors <- function(params) unlist(params$colors$timepoint)
diet_time_colors <- function(params) unlist(params$colors$diet_time)

format_p_value <- function(p) {
  fifelse(
    is.na(p),
    "p = NA",
    fifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  )
}

format_p_stars <- function(p) {
  fifelse(
    is.na(p), "ns",
    fifelse(p < 0.001, "***", fifelse(p < 0.01, "**", fifelse(p < 0.05, "*", "ns")))
  )
}

format_permanova_label <- function(dt, term_name, prefix) {
  row <- as.data.table(dt)[term == term_name]
  if (nrow(row) == 0) return(paste0(prefix, " p = NA"))
  r2_col <- intersect(c("R2", "R.squared"), names(row))[1]
  p_col <- intersect(c("Pr(>F)", "Pr..F."), names(row))[1]
  r2 <- if (!is.na(r2_col)) row[[r2_col]][1] else NA_real_
  p <- if (!is.na(p_col)) row[[p_col]][1] else NA_real_
  paste0(prefix, " R2 = ", sprintf("%.3f", r2), ", ", format_p_value(p))
}

save_plot <- function(plot, file_base, width = 7, height = 5) {
  ensure_dir(dirname(file_base))
  ggplot2::ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(paste0(file_base, ".svg"), plot, width = width, height = height, device = grDevices::svg)
}

clr_transform <- function(mat, pseudocount = 0.5) {
  mat <- as.matrix(mat)
  mat <- mat + pseudocount
  log_mat <- log(mat)
  sweep(log_mat, 1, rowMeans(log_mat), "-")
}

require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing R packages: ", paste(missing, collapse = ", "), ". Install downstream dependencies first.")
  }
}

sample_data_dt <- function(ps) {
  df <- as(phyloseq::sample_data(ps), "data.frame")
  dt <- as.data.table(df)
  if (!"sample_id" %in% names(dt)) dt[, sample_id := rownames(df)]
  setcolorder(dt, "sample_id")
  dt
}

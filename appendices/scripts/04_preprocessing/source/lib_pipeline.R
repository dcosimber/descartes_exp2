suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

project_root <- function() {
  normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "scripts"), ".."), mustWork = FALSE)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

script_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[startsWith(args, file_arg)][1])
  if (is.na(script_path)) {
    normalizePath(getwd())
  } else {
    normalizePath(file.path(dirname(script_path), ".."))
  }
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

load_params <- function(project = script_project_root()) {
  yaml::read_yaml(file.path(project, "config", "preprocessing_params.yaml"))
}

merge_lists <- function(base, override) {
  if (is.null(override)) return(base)
  for (nm in names(override)) {
    if (is.list(base[[nm]]) && is.list(override[[nm]]) &&
        !is.data.frame(base[[nm]]) && !is.data.frame(override[[nm]])) {
      base[[nm]] <- merge_lists(base[[nm]], override[[nm]])
    } else {
      base[[nm]] <- override[[nm]]
    }
  }
  base
}

load_run_params <- function(project = script_project_root(), run_id = NULL) {
  params <- load_params(project)
  if (!is.null(run_id) && !is.na(run_id) && nzchar(run_id)) {
    override_file <- file.path(project, "config", "run_overrides", paste0(run_id, ".yaml"))
    if (file.exists(override_file)) {
      params <- merge_lists(params, yaml::read_yaml(override_file))
    }
  }
  params
}

path_from_project <- function(project, path) {
  normalizePath(file.path(project, path), mustWork = FALSE)
}

path_from_descartes <- function(project, path) {
  normalizePath(file.path(project, path), mustWork = FALSE)
}

ensure_dir <- function(...) {
  dirs <- unlist(list(...), use.names = FALSE)
  for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

write_tsv <- function(x, file) {
  ensure_dir(dirname(file))
  data.table::fwrite(as.data.table(x), file, sep = "\t", na = "NA", quote = FALSE)
}

read_tsv <- function(file) {
  data.table::fread(file, sep = "\t", na.strings = c("NA", ""))
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_message <- function(...) {
  cat(sprintf("[%s] %s\n", timestamp(), paste(..., collapse = " ")))
}

split_subset <- function(x) {
  if (is.null(x) || isTRUE(x) || is.na(x) || x == "") character()
  else trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

sample_info_from_id <- function(sample_id) {
  if (grepl("^M[0-9]+$", sample_id)) {
    return(list(
      sample_type = "Mock", include_qc = "yes", include_biological_analysis = "no",
      diet_code = NA_character_, diet_label = "mock_community",
      hydrolysate_percent = NA_integer_, fishmeal_percent = NA_integer_,
      timepoint_day = NA_integer_, tank = NA_integer_, fish = NA_integer_,
      group = "Mock"
    ))
  }
  if (grepl("^[XYZ]3$", sample_id)) {
    code <- substr(sample_id, 1, 1)
    labels <- c(
      X = "commercial_control_feed",
      Y = "fishmeal_free_feed",
      Z = "blue_whiting_hydrolysate_30_feed"
    )
    return(list(
      sample_type = "Feed", include_qc = "yes", include_biological_analysis = "no",
      diet_code = code, diet_label = labels[[code]],
      hydrolysate_percent = ifelse(code == "Z", 30L, ifelse(code == "Y", 0L, 0L)),
      fishmeal_percent = ifelse(code == "X", 30L, 0L),
      timepoint_day = NA_integer_, tank = 3L, fish = NA_integer_,
      group = paste0(code, ".Feed")
    ))
  }
  m <- regexec("^([XYZ])([0-9]+)\\.([0-9]+)t([0-9]+)$", sample_id)
  parts <- regmatches(sample_id, m)[[1]]
  if (length(parts) == 5) {
    code <- parts[[2]]
    tank <- as.integer(parts[[3]])
    fish <- as.integer(parts[[4]])
    timepoint <- as.integer(parts[[5]])
    labels <- c(
      X = "commercial_control",
      Y = "fishmeal_free",
      Z = "blue_whiting_hydrolysate_30"
    )
    return(list(
      sample_type = "Intestine", include_qc = "yes", include_biological_analysis = "yes",
      diet_code = code, diet_label = labels[[code]],
      hydrolysate_percent = ifelse(code == "Z", 30L, 0L),
      fishmeal_percent = ifelse(code == "X", 30L, 0L),
      timepoint_day = timepoint, tank = tank, fish = fish,
      group = paste0(code, ".t", timepoint)
    ))
  }
  list(
    sample_type = "Unknown", include_qc = "yes", include_biological_analysis = "no",
    diet_code = NA_character_, diet_label = NA_character_,
    hydrolysate_percent = NA_integer_, fishmeal_percent = NA_integer_,
    timepoint_day = NA_integer_, tank = NA_integer_, fish = NA_integer_,
    group = "Unknown"
  )
}

expected_sample_ids <- function() {
  ids <- c("M1", "M2", "X3", "Y3", "Z3")
  for (diet in c("X", "Y", "Z")) {
    for (tank in 1:3) {
      for (fish in 1:5) {
        for (timepoint in c(7, 30, 90)) {
          ids <- c(ids, sprintf("%s%d.%dt%d", diet, tank, fish, timepoint))
        }
      }
    }
  }
  ids
}

read_run_metadata <- function(project, run_id) {
  file <- file.path(project, "metadata", paste0("sample_metadata_", run_id, ".tsv"))
  if (!file.exists(file)) file <- file.path(project, "metadata", "sample_metadata.tsv")
  read_tsv(file)
}

run_dir <- function(project, step, run_id) {
  file.path(project, "results", step, run_id)
}

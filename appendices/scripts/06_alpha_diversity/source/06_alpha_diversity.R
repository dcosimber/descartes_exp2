#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(phyloseq)
  library(ggplot2)
})

source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]])), "lib_downstream.R"))

args <- parse_cli(list(run_id = NULL, taxonomy_database = NULL))
project <- script_project_root()
params <- load_params(project)
run_id <- args$run_id %||% params$project$final_preprocessing_run
tax_db <- args$taxonomy_database %||% params$project$taxonomy_database

tables_dir <- stage_tables_dir(project, "06_alpha_diversity")
fig_dir <- stage_figures_dir(project, "06_alpha_diversity")
ensure_dir(tables_dir, fig_dir)

test_one <- function(dt, group_col, metric_col, subset_expr = TRUE) {
  d <- dt[eval(substitute(subset_expr), dt, parent.frame())]
  d <- d[!is.na(get(group_col))]
  groups <- unique(d[[group_col]])
  if (length(groups) < 2) return(data.table())
  global_test <- if (length(groups) == 2) {
    wt <- wilcox.test(d[[metric_col]] ~ d[[group_col]], exact = FALSE)
    data.table(test = "wilcoxon_rank_sum", group_variable = group_col, metric = metric_col, group1 = as.character(groups[1]), group2 = as.character(groups[2]), p_adj = wt$p.value)
  } else {
    kw <- kruskal.test(d[[metric_col]] ~ d[[group_col]])
    data.table(test = "kruskal_wallis", group_variable = group_col, metric = metric_col, group1 = "global", group2 = NA_character_, p_adj = kw$p.value)
  }
  if (length(groups) < 3) return(global_test)
  pairs <- pairwise.wilcox.test(d[[metric_col]], d[[group_col]], p.adjust.method = "BH", exact = FALSE)
  pair_dt <- as.data.table(as.table(pairs$p.value))
  setnames(pair_dt, c("group1", "group2", "p_adj"))
  pair_dt <- pair_dt[!is.na(p_adj)]
  pair_dt[, `:=`(test = "pairwise_wilcoxon", group_variable = group_col, metric = metric_col)]
  rbind(global_test, pair_dt, fill = TRUE)
}

compute_alpha_stats <- function(alpha, metric_cols) {
  stats <- list()
  for (metric in metric_cols) {
    stats[[paste(metric, "diet")]] <- test_one(alpha, "diet_code", metric)
    stats[[paste(metric, "hydro")]] <- test_one(alpha, "hydro_group", metric)
    stats[[paste(metric, "time")]] <- test_one(alpha, "timepoint", metric)
    for (diet in na.omit(unique(alpha$diet_code))) {
      tmp <- alpha[diet_code == diet]
      st <- test_one(tmp, "timepoint", metric)
      if (nrow(st)) st[, stratum := as.character(diet)]
      stats[[paste(metric, "time", diet)]] <- st
    }
    for (hydro in na.omit(unique(alpha$hydro_group))) {
      tmp <- alpha[hydro_group == hydro]
      st <- test_one(tmp, "timepoint", metric)
      if (nrow(st)) st[, stratum := as.character(hydro)]
      stats[[paste(metric, "time", hydro)]] <- st
    }
    for (tp in na.omit(unique(alpha$timepoint))) {
      tmp <- alpha[timepoint == tp]
      st <- test_one(tmp, "diet_code", metric)
      if (nrow(st)) st[, stratum := as.character(tp)]
      stats[[paste(metric, "diet", tp)]] <- st
      st <- test_one(tmp, "hydro_group", metric)
      if (nrow(st)) st[, stratum := as.character(tp)]
      stats[[paste(metric, "hydro", tp)]] <- st
    }
  }
  rbindlist(stats, fill = TRUE)
}

alpha_metric_label <- function(metric) {
  labels <- c(Observed = "Observed richness", Shannon = "Shannon index", Simpson = "Simpson index", Pielou = "Pielou evenness")
  metric <- as.character(metric)
  out <- unname(labels[metric])
  out[is.na(out)] <- metric[is.na(out)]
  out
}

alpha_level_label <- function(tax_level) {
  labels <- c(ASV = "ASV", Genus = "Genus")
  labels[[tax_level]] %||% tax_level
}

darken_colors <- function(cols, factor = 0.58) {
  rgb_mat <- grDevices::col2rgb(cols)
  dark_rgb <- matrix(pmax(0, round(rgb_mat * factor)), nrow = 3)
  out <- grDevices::rgb(dark_rgb[1, ], dark_rgb[2, ], dark_rgb[3, ], maxColorValue = 255)
  names(out) <- names(cols)
  out
}

alpha_bracket_data <- function(plot_dt, stats_dt, group_var, facet_col, x_levels) {
  labels <- stats_dt[get("group_variable") == group_var & !is.na(stratum)]
  labels <- labels[stratum %in% as.character(plot_dt[[facet_col]])]
  if (nrow(labels) == 0) return(data.table())
  labels[, (facet_col) := factor(stratum, levels = levels(plot_dt[[facet_col]]))]
  labels[, `:=`(
    xstart = fifelse(test == "kruskal_wallis", 1, match(group1, x_levels)),
    xend = fifelse(test == "kruskal_wallis", length(x_levels), match(group2, x_levels))
  )]
  labels <- labels[!is.na(xstart) & !is.na(xend)]
  labels[, `:=`(
    xleft = pmin(xstart, xend),
    xright = pmax(xstart, xend)
  )]
  labels[, bracket_type := fifelse(test == "kruskal_wallis", 2L, 1L)]
  setorderv(labels, c("metric", facet_col, "bracket_type", "xleft", "xright"))
  labels[, bracket_rank := seq_len(.N), by = c("metric", facet_col)]
  ypos <- plot_dt[, .(
    ymin = min(value, na.rm = TRUE),
    ymax = max(value, na.rm = TRUE)
  ), by = c("metric", facet_col)]
  ypos[, yrange := pmax(ymax - ymin, ymax * 0.08, 0.1)]
  labels <- merge(labels, ypos, by = c("metric", facet_col))
  labels[, `:=`(
    y = ymax + yrange * (0.11 + 0.12 * bracket_rank),
    tick = yrange * 0.025,
    label = format_p_stars(p_adj)
  )]
  labels[]
}

alpha_nested_bracket_data <- function(plot_dt, stats_dt, x_lookup) {
  global <- stats_dt[
    group_variable == "hydro_group" &
      test == "wilcoxon_rank_sum" &
      !is.na(stratum) &
      metric %in% unique(plot_dt$metric) &
      stratum %in% as.character(unique(plot_dt$timepoint))
  ]
  if (nrow(global) > 0) {
    global[, `:=`(
      timepoint = factor(stratum, levels = levels(plot_dt$timepoint)),
      group1 = "Ctrl",
      group2 = "Hydrolysate",
      xleft = unname(x_lookup["Ctrl"]),
      xright = unname(x_lookup["BL30"]),
      bracket_type = 2L
    )]
  }
  pairwise <- stats_dt[
    group_variable == "diet_code" &
      test == "pairwise_wilcoxon" &
      !is.na(stratum) &
      metric %in% unique(plot_dt$metric) &
      stratum %in% as.character(unique(plot_dt$timepoint))
  ]
  pairwise <- pairwise[
    (as.character(group1) == "BL15" & as.character(group2) == "BL30") |
      (as.character(group1) == "BL30" & as.character(group2) == "BL15")
  ]
  if (nrow(pairwise) > 0) {
    pairwise[, `:=`(
      timepoint = factor(stratum, levels = levels(plot_dt$timepoint)),
      xstart = unname(x_lookup[as.character(group1)]),
      xend = unname(x_lookup[as.character(group2)]),
      bracket_type = 1L
    )]
    pairwise <- pairwise[!is.na(xstart) & !is.na(xend)]
    pairwise[, `:=`(
      xleft = pmin(xstart, xend),
      xright = pmax(xstart, xend)
    )]
  }
  labels <- rbindlist(list(pairwise, global), fill = TRUE)
  if (nrow(labels) == 0) return(data.table())
  labels <- labels[!is.na(xleft) & !is.na(xright)]
  setorderv(labels, c("metric", "timepoint", "bracket_type", "xleft", "xright"))
  ypos <- plot_dt[, .(
    ymin = min(value, na.rm = TRUE),
    ymax = max(value, na.rm = TRUE)
  ), by = .(metric, timepoint)]
  ypos[, yrange := pmax(ymax - ymin, ymax * 0.08, 0.1)]
  labels <- merge(labels, ypos, by = c("metric", "timepoint"))
  labels[, `:=`(
    y = ymax + fifelse(bracket_type == 1L, yrange * 0.10, yrange * 0.30),
    tick = yrange * 0.018,
    label = format_p_stars(p_adj)
  )]
  labels[]
}

add_alpha_brackets <- function(plot, bracket_dt) {
  if (nrow(bracket_dt) == 0) return(plot)
  plot +
    geom_segment(
      data = bracket_dt,
      aes(x = xleft, xend = xright, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.28,
      color = "grey15"
    ) +
    geom_segment(
      data = bracket_dt,
      aes(x = xleft, xend = xleft, y = y, yend = y - tick),
      inherit.aes = FALSE,
      linewidth = 0.28,
      color = "grey15"
    ) +
    geom_segment(
      data = bracket_dt,
      aes(x = xright, xend = xright, y = y, yend = y - tick),
      inherit.aes = FALSE,
      linewidth = 0.28,
      color = "grey15"
    ) +
    geom_text(
      data = bracket_dt,
      aes(x = (xleft + xright) / 2, y = y + tick * 0.8, label = label),
      inherit.aes = FALSE,
      size = 2.05,
      color = "grey15",
      vjust = 0
    )
}

save_alpha_plot <- function(plot, base_name, suffix, width, height, legacy = FALSE) {
  comparison_dir <- sub("^alpha_by_", "", base_name)
  save_plot(plot, file.path(fig_dir, comparison_dir, paste0(base_name, "_", suffix)), width = width, height = height)
  if (legacy) save_plot(plot, file.path(fig_dir, comparison_dir, base_name), width = width, height = height)
}

make_diet_time_plot <- function(alpha_long, stats_dt, tax_level) {
  nested_x <- c(Ctrl = 1, BL15 = 1.82, BL30 = 2.18)
  plot_dt <- copy(alpha_long)
  plot_dt[, x_plot := unname(nested_x[as.character(diet_code)])]
  plot_dt <- plot_dt[!is.na(x_plot)]
  plot_dt[, diet_code := factor(diet_code, levels = params$design$diet_levels)]
  brackets <- alpha_nested_bracket_data(plot_dt, stats_dt, nested_x)

  p <- ggplot(plot_dt, aes(x_plot, value, fill = diet_code, group = diet_code)) +
    geom_boxplot(width = 0.28, outlier.shape = NA, alpha = 0.72, linewidth = 0.35, color = "black") +
    geom_jitter(
      aes(color = diet_code),
      width = 0.035,
      size = 1.15,
      alpha = 0.48,
      shape = 21,
      stroke = 0.32,
      show.legend = FALSE
    ) +
    facet_grid(metric ~ timepoint, scales = "free_y", labeller = labeller(metric = alpha_metric_label), switch = "y") +
    scale_fill_manual(values = diet_colors(params), drop = FALSE) +
    scale_color_manual(values = darken_colors(diet_colors(params)), drop = FALSE) +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Ctrl", "Hydrolysate"),
      limits = c(0.65, 2.45),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.32))) +
    labs(x = NULL, y = NULL, fill = "Diet", color = "Diet", title = alpha_level_label(tax_level)) +
    theme_publication(base_size = 8) +
    coord_cartesian(clip = "off") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 9),
      strip.text = element_text(size = 8.5, color = "grey10"),
      strip.text.y.left = element_text(angle = 90, size = 8.5, face = "bold"),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      axis.text.x = element_text(color = "grey20", size = 7.5),
      axis.text.y = element_text(color = "grey20", size = 7.5),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 8),
      legend.text = element_text(size = 7.5),
      legend.key.size = grid::unit(0.35, "cm"),
      panel.border = element_rect(color = "grey25", fill = NA, linewidth = 0.35),
      panel.spacing = grid::unit(0.45, "lines"),
      plot.margin = margin(5, 8, 5, 5)
    )
  add_alpha_brackets(p, brackets)
}

make_hydro_time_plot <- function(alpha_long, stats_dt, tax_level) {
  plot_dt <- alpha_long[!is.na(hydro_group)]
  brackets <- alpha_bracket_data(plot_dt, stats_dt, "hydro_group", "timepoint", params$design$hydro_group_levels)
  p <- ggplot(plot_dt, aes(hydro_group, value, fill = hydro_group, group = hydro_group)) +
    geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.72, color = "black") +
    geom_jitter(aes(color = hydro_group), width = 0.10, size = 1.2, alpha = 0.48, shape = 21, stroke = 0.32, show.legend = FALSE) +
    facet_grid(metric ~ timepoint, scales = "free_y", labeller = labeller(metric = alpha_metric_label), switch = "y") +
    scale_fill_manual(values = hydro_colors(params), drop = FALSE) +
    scale_color_manual(values = darken_colors(hydro_colors(params)), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.32))) +
    labs(x = NULL, y = NULL, fill = "Group", color = "Group", title = paste0(alpha_level_label(tax_level), ": control vs hydrolysate")) +
    theme_publication(base_size = 8) +
    coord_cartesian(clip = "off") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 9),
      strip.text.y.left = element_text(angle = 90, size = 8.5, face = "bold"),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      legend.position = "bottom",
      plot.margin = margin(6, 12, 6, 6)
    )
  add_alpha_brackets(p, brackets)
}

make_time_diet_plot <- function(alpha_long, stats_dt, tax_level) {
  time_levels <- sprintf("D%02d", params$design$timepoint_levels)
  plot_dt <- copy(alpha_long)
  plot_dt[, timepoint := factor(timepoint, levels = time_levels)]
  plot_dt[, diet_code := factor(diet_code, levels = params$design$diet_levels)]
  brackets <- alpha_bracket_data(plot_dt, stats_dt, "timepoint", "diet_code", time_levels)

  p <- ggplot(plot_dt, aes(timepoint, value, fill = timepoint, group = timepoint)) +
    geom_boxplot(width = 0.60, outlier.shape = NA, alpha = 0.72, color = "black") +
    geom_jitter(aes(color = timepoint), width = 0.10, size = 1.15, alpha = 0.48, shape = 21, stroke = 0.32, show.legend = FALSE) +
    facet_grid(metric ~ diet_code, scales = "free_y", labeller = labeller(metric = alpha_metric_label), switch = "y") +
    scale_fill_manual(values = time_colors(params), drop = FALSE) +
    scale_color_manual(values = darken_colors(time_colors(params)), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.46))) +
    labs(x = "Timepoint", y = NULL, fill = "Timepoint", color = "Timepoint", title = alpha_level_label(tax_level)) +
    theme_publication(base_size = 8) +
    coord_cartesian(clip = "off") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 9),
      strip.text = element_text(size = 8.5, color = "grey10"),
      strip.text.y.left = element_text(angle = 90, size = 8.5, face = "bold"),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      axis.text.x = element_text(color = "grey20", size = 7.5),
      axis.text.y = element_text(color = "grey20", size = 7.5),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 8),
      legend.text = element_text(size = 7.5),
      legend.key.size = grid::unit(0.35, "cm"),
      panel.border = element_rect(color = "grey25", fill = NA, linewidth = 0.35),
      panel.spacing = grid::unit(0.45, "lines"),
      plot.margin = margin(5, 8, 5, 5)
    )
  add_alpha_brackets(p, brackets)
}

feed_colors <- function(params) {
  c(Ctrl = unname(diet_colors(params)["Ctrl"]), Hydrolysate = unname(hydro_colors(params)["Hydrolysate"]), Feed = "#7A7A7A")
}

alpha_feed_bracket_data <- function(plot_dt, stats_dt) {
  labels <- stats_dt[group_variable == "ctrl_hydro_feed_group"]
  if (nrow(labels) == 0) return(data.table())
  x_levels <- levels(plot_dt$ctrl_hydro_feed_group)
  labels[, `:=`(
    xstart = fifelse(test == "kruskal_wallis", 1, match(group1, x_levels)),
    xend = fifelse(test == "kruskal_wallis", length(x_levels), match(group2, x_levels))
  )]
  labels <- labels[!is.na(xstart) & !is.na(xend)]
  labels[, `:=`(
    xleft = pmin(xstart, xend),
    xright = pmax(xstart, xend),
    bracket_type = fifelse(test == "kruskal_wallis", 2L, 1L)
  )]
  setorderv(labels, c("tax_level", "metric", "bracket_type", "xleft", "xright"))
  labels[, bracket_rank := seq_len(.N), by = .(tax_level, metric)]
  ypos <- plot_dt[, .(
    ymin = min(value, na.rm = TRUE),
    ymax = max(value, na.rm = TRUE)
  ), by = .(tax_level, metric)]
  ypos[, yrange := pmax(ymax - ymin, ymax * 0.08, 0.1)]
  labels <- merge(labels, ypos, by = c("tax_level", "metric"))
  labels[, `:=`(
    y = ymax + yrange * (0.10 + 0.105 * bracket_rank),
    tick = yrange * 0.020,
    label = format_p_stars(p_adj)
  )]
  labels[]
}

make_ctrl_hydro_feed_plot <- function(alpha_feed_long, feed_stats) {
  plot_dt <- copy(alpha_feed_long)
  plot_dt[, ctrl_hydro_feed_group := factor(ctrl_hydro_feed_group, levels = c("Ctrl", "Hydrolysate", "Feed"))]
  brackets <- alpha_feed_bracket_data(plot_dt, feed_stats)
  p <- ggplot(plot_dt, aes(ctrl_hydro_feed_group, value, fill = ctrl_hydro_feed_group, group = ctrl_hydro_feed_group)) +
    geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.72, color = "black") +
    geom_jitter(aes(color = ctrl_hydro_feed_group), width = 0.09, size = 1.15, alpha = 0.48, shape = 21, stroke = 0.32, show.legend = FALSE) +
    facet_grid(metric ~ tax_level, scales = "free_y", labeller = labeller(metric = alpha_metric_label), switch = "y") +
    scale_fill_manual(values = feed_colors(params), drop = FALSE) +
    scale_color_manual(values = darken_colors(feed_colors(params)), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.42))) +
    labs(x = NULL, y = NULL, fill = "Group", color = "Group") +
    theme_publication(base_size = 8) +
    coord_cartesian(clip = "off") +
    theme(
      strip.text = element_text(size = 8.5, color = "grey10"),
      strip.text.y.left = element_text(angle = 90, size = 8.5, face = "bold"),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      axis.text.x = element_text(color = "grey20", size = 7.3),
      axis.text.y = element_text(color = "grey20", size = 7.5),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 8),
      legend.text = element_text(size = 7.5),
      legend.key.size = grid::unit(0.35, "cm"),
      panel.border = element_rect(color = "grey25", fill = NA, linewidth = 0.35),
      panel.spacing = grid::unit(0.45, "lines"),
      plot.margin = margin(5, 8, 5, 5)
    )
  add_alpha_brackets(p, brackets)
}

alpha_objects <- list(ASV = "ps_biological_final", Genus = "ps_biological_genus")
all_alpha <- list()
all_stats <- list()
all_alpha_feed <- list()
all_feed_stats <- list()
ps_final_all <- read_object(project, "02_filtering", "ps_final", tax_db, run_id)

for (tax_level in names(alpha_objects)) {
  ps <- read_object(project, "02_filtering", alpha_objects[[tax_level]], tax_db, run_id)
  suffix <- tolower(tax_level)
  legacy <- identical(tax_level, "ASV")

  # Alpha diversity is estimated on the filtered downstream objects; phyloseq warns
  # about missing singletons, which is expected after prevalence/abundance filtering.
  alpha_measures_base <- setdiff(params$analysis$alpha_measures, "Pielou")
  alpha <- as.data.table(suppressWarnings(estimate_richness(ps, measures = alpha_measures_base)), keep.rownames = "sample_id")
  if ("Pielou" %in% params$analysis$alpha_measures) {
    alpha[, Pielou := fifelse(Observed > 1 & Shannon > 0, Shannon / log(Observed), 0)]
  }
  alpha <- merge(alpha, sample_data_dt(ps), by = "sample_id")
  alpha[, tax_level := tax_level]
  write_tsv(alpha, file.path(tables_dir, paste0("alpha_diversity_", suffix, ".tsv")))

  metric_cols <- intersect(params$analysis$alpha_measures, names(alpha))
  alpha_long <- melt(alpha, id.vars = setdiff(names(alpha), metric_cols), measure.vars = metric_cols, variable.name = "metric", value.name = "value")
  alpha_long[, metric := factor(metric, levels = params$analysis$alpha_measures)]
  alpha_long[, timepoint := factor(timepoint, levels = sprintf("D%02d", params$design$timepoint_levels))]
  alpha_long[, diet_code := factor(diet_code, levels = params$design$diet_levels)]
  alpha_long[, hydro_group := factor(hydro_group, levels = params$design$hydro_group_levels)]

  stats_dt <- compute_alpha_stats(alpha, metric_cols)
  stats_dt[, tax_level := tax_level]
  write_tsv(stats_dt, file.path(tables_dir, paste0("alpha_tests_", suffix, ".tsv")))

  p_diet <- make_diet_time_plot(alpha_long, stats_dt, tax_level)
  save_alpha_plot(p_diet, "alpha_by_diet_time", suffix, width = 6.35, height = 8.25, legacy = legacy)

  p_hydro <- make_hydro_time_plot(alpha_long, stats_dt, tax_level)
  save_alpha_plot(p_hydro, "alpha_by_hydro_time", suffix, width = 6.0, height = 8.25, legacy = legacy)

  p_time <- make_time_diet_plot(alpha_long, stats_dt, tax_level)
  save_alpha_plot(p_time, "alpha_by_time_diet", suffix, width = 6.45, height = 8.25, legacy = legacy)

  feed_ps <- subset_samples(ps_final_all, (sample_type == params$design$biological_sample_type & include_biological_analysis == "yes") | sample_type == "Feed")
  feed_ps <- prune_taxa(taxa_sums(feed_ps) > 0, feed_ps)
  if (identical(tax_level, "Genus")) {
    feed_ps <- tax_glom(feed_ps, taxrank = "Genus", NArm = FALSE)
    feed_ps <- prune_taxa(taxa_sums(feed_ps) > 0, feed_ps)
  }
  alpha_feed_measures_base <- setdiff(params$analysis$alpha_measures, "Pielou")
  alpha_feed <- as.data.table(suppressWarnings(estimate_richness(feed_ps, measures = alpha_feed_measures_base)), keep.rownames = "sample_id")
  if ("Pielou" %in% params$analysis$alpha_measures) {
    alpha_feed[, Pielou := fifelse(Observed > 1 & Shannon > 0, Shannon / log(Observed), 0)]
  }
  alpha_feed <- merge(alpha_feed, sample_data_dt(feed_ps), by = "sample_id")
  alpha_feed[, tax_level := tax_level]
  alpha_feed[, ctrl_hydro_feed_group := fifelse(sample_type == "Feed", "Feed", as.character(hydro_group))]
  alpha_feed[, ctrl_hydro_feed_group := factor(ctrl_hydro_feed_group, levels = c("Ctrl", "Hydrolysate", "Feed"))]

  feed_metric_cols <- intersect(params$analysis$alpha_measures, names(alpha_feed))
  feed_stats <- rbindlist(lapply(feed_metric_cols, function(metric) test_one(alpha_feed, "ctrl_hydro_feed_group", metric)), fill = TRUE)
  feed_stats[, tax_level := tax_level]
  alpha_feed_long <- melt(alpha_feed, id.vars = setdiff(names(alpha_feed), feed_metric_cols), measure.vars = feed_metric_cols, variable.name = "metric", value.name = "value")
  alpha_feed_long[, metric := factor(metric, levels = params$analysis$alpha_measures)]
  all_alpha_feed[[tax_level]] <- alpha_feed
  all_feed_stats[[tax_level]] <- feed_stats

  all_alpha[[tax_level]] <- alpha
  all_stats[[tax_level]] <- stats_dt
}

alpha_all <- rbindlist(all_alpha, fill = TRUE)
stats_all <- rbindlist(all_stats, fill = TRUE)
write_tsv(alpha_all, file.path(tables_dir, "alpha_diversity.tsv"))
write_tsv(stats_all, file.path(tables_dir, "alpha_tests.tsv"))

alpha_feed_all <- rbindlist(all_alpha_feed, fill = TRUE)
feed_stats_all <- rbindlist(all_feed_stats, fill = TRUE)
write_tsv(alpha_feed_all, file.path(tables_dir, "alpha_diversity_ctrl_hydro_feed.tsv"))
write_tsv(feed_stats_all, file.path(tables_dir, "alpha_tests_ctrl_hydro_feed.tsv"))

alpha_feed_long_all <- melt(alpha_feed_all, id.vars = setdiff(names(alpha_feed_all), params$analysis$alpha_measures), measure.vars = params$analysis$alpha_measures, variable.name = "metric", value.name = "value")
alpha_feed_long_all[, metric := factor(metric, levels = params$analysis$alpha_measures)]
p_feed <- make_ctrl_hydro_feed_plot(alpha_feed_long_all, feed_stats_all)
save_plot(p_feed, file.path(fig_dir, "ctrl_hydro_feed", "alpha_ctrl_hydro_feed"), width = 5.2, height = 7.8)

write_xlsx(
  list(alpha_values = alpha_all, alpha_tests = stats_all, alpha_ctrl_hydro_feed = alpha_feed_all, tests_ctrl_hydro_feed = feed_stats_all),
  file.path(tables_dir, "alpha_diversity.xlsx")
)

log_message("Alpha diversity complete:", run_id)

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dada2)
})

source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]])), "lib_pipeline.R"))

args <- parse_cli(list(run_id = "run_01_default"))
project <- script_project_root()
run_id <- args$run_id
params <- load_run_params(project, run_id)
dada_params <- params$dada2$dada
merge_params <- params$dada2$merge_pairs
chimera_params <- params$dada2$chimera

filter_stats <- read_tsv(file.path(run_dir(project, "02_filter_trim", run_id), "dada2_filter_stats.tsv"))
err_f <- readRDS(file.path(run_dir(project, "03_error_models", run_id), "errF.rds"))
err_r <- readRDS(file.path(run_dir(project, "03_error_models", run_id), "errR.rds"))

outdir <- run_dir(project, "04_denoise_merge", run_id)
chimera_dir <- run_dir(project, "05_chimera", run_id)
ensure_dir(outdir, chimera_dir, file.path(project, "logs", "commands"))

cmd_log <- file.path(project, "logs", "commands", paste0("06_denoise_merge_chimera_", run_id, ".commands.R"))
writeLines(c(capture.output(str(dada_params)), capture.output(str(merge_params)), capture.output(str(chimera_params))), cmd_log)

sample_names <- filter_stats$sample_id
log_message("Dereplicating filtered reads.")
derep_f <- derepFastq(filter_stats$filt_f, verbose = TRUE)
derep_r <- derepFastq(filter_stats$filt_r, verbose = TRUE)
names(derep_f) <- sample_names
names(derep_r) <- sample_names

log_message("Denoising forward reads.")
dada_f <- dada(
  derep_f,
  err = err_f,
  multithread = TRUE,
  pool = isTRUE(dada_params$pool),
  selfConsist = isTRUE(dada_params$self_consist)
)
log_message("Denoising reverse reads.")
dada_r <- dada(
  derep_r,
  err = err_r,
  multithread = TRUE,
  pool = isTRUE(dada_params$pool),
  selfConsist = isTRUE(dada_params$self_consist)
)

log_message("Merging paired reads.")
mergers <- mergePairs(
  dadaF = dada_f,
  derepF = derep_f,
  dadaR = dada_r,
  derepR = derep_r,
  minOverlap = as.integer(merge_params$min_overlap),
  maxMismatch = as.integer(merge_params$max_mismatch),
  verbose = TRUE
)

seqtab <- makeSequenceTable(mergers)
log_message("Removing chimeras.")
seqtab_nochim <- removeBimeraDenovo(
  seqtab,
  method = chimera_params$method,
  multithread = TRUE,
  verbose = TRUE
)

get_n <- function(x) sum(getUniques(x))
track <- data.table(
  sample_id = sample_names,
  filtered_input = filter_stats$reads.in,
  filtered_output = filter_stats$reads.out,
  denoised_forward = sapply(dada_f, get_n),
  denoised_reverse = sapply(dada_r, get_n),
  merged = rowSums(seqtab),
  nonchim = rowSums(seqtab_nochim)
)
track[, `:=`(
  filter_retention_pct = filtered_output / filtered_input * 100,
  merge_retention_pct = merged / filtered_input * 100,
  nonchim_retention_pct = nonchim / filtered_input * 100,
  chimera_loss_pct = (merged - nonchim) / filtered_input * 100
)]

saveRDS(derep_f, file.path(outdir, "derepF.rds"))
saveRDS(derep_r, file.path(outdir, "derepR.rds"))
saveRDS(dada_f, file.path(outdir, "dadaF.rds"))
saveRDS(dada_r, file.path(outdir, "dadaR.rds"))
saveRDS(mergers, file.path(outdir, "mergers.rds"))
saveRDS(seqtab, file.path(outdir, "seqtab.rds"))
saveRDS(seqtab_nochim, file.path(chimera_dir, "seqtab_nochim.rds"))

write_tsv(track, file.path(chimera_dir, "dada2_merge_chimera_stats.tsv"))
write_tsv(track, file.path(project, "results", "08_qc_tables", paste0(run_id, "_dada2_merge_chimera_stats.tsv")))

seq_lengths <- data.table(sequence = colnames(seqtab_nochim), length = nchar(colnames(seqtab_nochim)), total_abundance = colSums(seqtab_nochim))
write_tsv(seq_lengths, file.path(chimera_dir, "asv_length_distribution.tsv"))
log_message("Denoise/merge/chimera complete. ASVs:", ncol(seqtab_nochim), "reads:", sum(seqtab_nochim))

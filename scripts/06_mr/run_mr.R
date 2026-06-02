#!/usr/bin/env Rscript
# 06_mr/run_mr.R
# Two-sample MR for all cohorts: Wald (n=1) or IVW-MRE (n>=2).
# BH-FDR correction within each cohort.
#
# Usage:
#   Rscript scripts/06_mr/run_mr.R [--cohort ARIC_EA] [--limit 100]

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
if (length(script_arg) == 0 || is.na(script_arg)) {
  repo_root <- normalizePath(".")
} else {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
  repo_root  <- normalizePath(file.path(script_dir, "..", ".."))
}

source(file.path(repo_root, "scripts", "rlib", "logging.R"))
source(file.path(repo_root, "scripts", "rlib", "checkpoint.R"))
source(file.path(repo_root, "scripts", "rlib", "progress.R"))
source(file.path(repo_root, "scripts", "rlib", "mr_methods.R"))
source(file.path(repo_root, "scripts", "rlib", "config.R"))

setup_logger("06_mr")

# ── CLI args ──────────────────────────────────────────────────────────────────
args     <- commandArgs(trailingOnly = TRUE)
cohorts_all <- pipeline_cohorts()
cohort_arg  <- "all"
limit_arg   <- Inf

i <- 1
while (i <= length(args)) {
  if (args[i] == "--cohort") { cohort_arg <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--limit") { limit_arg <- as.integer(args[i+1]); i <- i + 2 }
  else i <- i + 1
}
run_cohorts <- if (cohort_arg == "all") cohorts_all else cohort_arg

# ── Main ───────────────────────────────────────────────────────────────────────
run_cohort_mr <- function(cohort) {
  harm_dir  <- file.path("processed_data", cohort, "harmonised")
  out_path  <- file.path("processed_data", cohort, "mr_results.tsv")
  state_path <- file.path("processed_data", cohort, "_state_06.rds")

  tsv_files <- list.files(harm_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (is.finite(limit_arg)) tsv_files <- head(tsv_files, limit_arg)

  cp <- checkpoint_load(state_path)
  todo <- tsv_files[basename(tsv_files) %in%
                     paste0(remaining(tools::file_path_sans_ext(basename(tsv_files)), cp), ".tsv")]

  log_info("%s: %d files, %d remaining", cohort, length(tsv_files), length(todo))

  results <- list()
  tick <- serial_progress(length(todo), label = paste(cohort, "MR"))

  for (f in todo) {
    tick()
    seqid <- tools::file_path_sans_ext(basename(f))
    tryCatch({
      harm <- fread(f, sep = "\t")
      if (nrow(harm) == 0) { cp <- checkpoint_mark(cp, seqid); next }

      res <- run_protein_mr(harm)
      if (!is.null(res)) {
        results[[seqid]] <- data.frame(
          cohort  = cohort,
          seqid   = res$seqid,
          gene    = if ("gene" %in% names(harm)) harm$gene[1] else NA,
          n_snps  = res$n_snps,
          method  = res$method,
          beta    = res$beta,
          se      = res$se,
          pval    = res$pval,
          OR      = res$OR,
          OR_lo95 = res$OR_lo95,
          OR_hi95 = res$OR_hi95,
          stringsAsFactors = FALSE
        )
      }
      cp <- checkpoint_mark(cp, seqid)
    }, error = function(e) {
      log_warn("%s %s: MR failed — %s", cohort, seqid, conditionMessage(e))
      cp <<- checkpoint_mark_failed(cp, seqid, conditionMessage(e))
    })
  }

  if (length(results) == 0) {
    log_info("%s: no MR results", cohort)
    return(invisible(NULL))
  }

  dt <- rbindlist(results)
  # BH-FDR within this cohort
  dt[, fdr_q   := bh_fdr(pval)]
  dt[, fdr_pass := fdr_q < 0.05]

  # Append to existing results if resuming
  if (file.exists(out_path)) {
    existing <- fread(out_path, sep = "\t")
    dt <- rbindlist(list(existing[!seqid %in% dt$seqid], dt), fill = TRUE)
    # Recompute FDR over full set
    dt[, fdr_q   := bh_fdr(pval)]
    dt[, fdr_pass := fdr_q < 0.05]
  }

  fwrite(dt, out_path, sep = "\t")
  log_info("%s: wrote %d results (%d FDR<0.05) → %s",
           cohort, nrow(dt), sum(dt$fdr_pass, na.rm = TRUE), out_path)
  invisible(dt)
}

all_results <- lapply(run_cohorts, run_cohort_mr)

# ── Combined cross-cohort table ───────────────────────────────────────────────
cohort_files <- file.path("processed_data", cohorts_all, "mr_results.tsv")
cohort_files <- cohort_files[file.exists(cohort_files)]
if (length(cohort_files) > 0) {
  combined <- rbindlist(lapply(cohort_files, fread), fill = TRUE)
  combined[, fdr_pass_any := any(fdr_pass), by = seqid]
  fwrite(combined, "processed_data/mr_all_cohorts.tsv", sep = "\t")
  log_info("Combined: %d rows, %d proteins with FDR<0.05 in any cohort",
           nrow(combined), uniqueN(combined[fdr_pass_any == TRUE, seqid]))
}

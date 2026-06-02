#!/usr/bin/env Rscript
# 07_sensitivity/run_sensitivity.R
# Sensitivity analyses: Q/I², weighted median/mode, Egger, Steiger filtering.
# Emits `passes_sensitivity` column.
#
# Usage:
#   Rscript scripts/07_sensitivity/run_sensitivity.R [--cohort ARIC_EA] [--limit 100]

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
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

setup_logger("07_sensitivity")

# ── CLI args ──────────────────────────────────────────────────────────────────
args        <- commandArgs(trailingOnly = TRUE)
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
run_cohort_sensitivity <- function(cohort) {
  harm_dir   <- file.path("processed_data", cohort, "harmonised")
  mr_path    <- file.path("processed_data", cohort, "mr_results.tsv")
  out_path   <- file.path("processed_data", cohort, "sensitivity.tsv")
  state_path <- file.path("processed_data", cohort, "_state_07.rds")

  if (!file.exists(mr_path)) {
    log_warn("%s: mr_results.tsv not found — run 06_mr first", cohort)
    return(invisible(NULL))
  }

  mr_res <- fread(mr_path, sep = "\t")
  # Only compute sensitivity for proteins with ≥2 instruments
  candidates <- mr_res[n_snps >= 2, seqid]
  tsv_files  <- file.path(harm_dir, paste0(candidates, ".tsv"))
  tsv_files  <- tsv_files[file.exists(tsv_files)]
  if (is.finite(limit_arg)) tsv_files <- head(tsv_files, limit_arg)

  cp <- checkpoint_load(state_path)
  todo <- tsv_files[basename(tsv_files) %in%
                     paste0(remaining(tools::file_path_sans_ext(basename(tsv_files)), cp), ".tsv")]

  log_info("%s: %d candidates, %d remaining sensitivity analyses", cohort, length(tsv_files), length(todo))

  results <- list()
  tick <- serial_progress(length(todo), label = paste(cohort, "sensitivity"))

  for (f in todo) {
    tick()
    seqid <- tools::file_path_sans_ext(basename(f))
    tryCatch({
      harm <- fread(f, sep = "\t")
      if (nrow(harm) == 0) { cp <- checkpoint_mark(cp, seqid); next }

      sens <- run_sensitivity(harm)
      results[[seqid]] <- as.data.frame(sens)
      cp <- checkpoint_mark(cp, seqid)
    }, error = function(e) {
      log_warn("%s %s: sensitivity failed — %s", cohort, seqid, conditionMessage(e))
      cp <<- checkpoint_mark_failed(cp, seqid, conditionMessage(e))
    })
  }

  if (length(results) == 0) {
    log_info("%s: no sensitivity results", cohort)
    return(invisible(NULL))
  }

  dt <- rbindlist(results, fill = TRUE)

  # Append to existing if resuming
  if (file.exists(out_path)) {
    existing <- fread(out_path, sep = "\t")
    dt <- rbindlist(list(existing[!seqid %in% dt$seqid], dt), fill = TRUE)
  }

  fwrite(dt, out_path, sep = "\t")
  n_pass <- sum(dt$passes_sensitivity, na.rm = TRUE)
  log_info("%s: %d proteins, %d pass sensitivity → %s",
           cohort, nrow(dt), n_pass, out_path)
  invisible(dt)
}

lapply(run_cohorts, run_cohort_sensitivity)

# ── Combined cross-cohort sensitivity table + enrich mr_all_cohorts ──────────
sens_cohorts <- cohorts_all
sens_files  <- file.path("processed_data", sens_cohorts, "sensitivity.tsv")
exists_mask <- file.exists(sens_files)
sens_files  <- sens_files[exists_mask]
sens_cohorts <- sens_cohorts[exists_mask]

if (length(sens_files) > 0) {
  sens_all <- rbindlist(
    mapply(function(f, cohort) {
      dt <- fread(f, sep = "\t")
      dt[, cohort := cohort]
      dt
    }, sens_files, sens_cohorts, SIMPLIFY = FALSE),
    fill = TRUE
  )
  fwrite(sens_all, "processed_data/sensitivity_all_cohorts.tsv", sep = "\t")
  log_info("Combined sensitivity: %d rows across %d cohorts",
           nrow(sens_all), length(sens_files))

  # Enrich mr_all_cohorts.tsv: left-join sensitivity columns by (cohort, seqid)
  mr_path <- "processed_data/mr_all_cohorts.tsv"
  if (file.exists(mr_path)) {
    mr_all <- fread(mr_path, sep = "\t")
    check <- merge(
      mr_all[, .(cohort, seqid, n_snps_mr = n_snps)],
      sens_all[, .(cohort, seqid, n_snps_sens = n_snps)],
      by = c("cohort", "seqid")
    )
    bad <- check[n_snps_mr != n_snps_sens]
    if (nrow(bad) > 0) {
      stop(sprintf(
        "n_snps mismatch between MR and sensitivity for %d proteins (e.g. %s/%s: MR=%d, sens=%d)",
        nrow(bad), bad$cohort[1], bad$seqid[1], bad$n_snps_mr[1], bad$n_snps_sens[1]
      ))
    }

    # n_snps is already in mr_all; drop it from sens to avoid a conflict.
    join_cols <- setdiff(names(sens_all), "n_snps")
    # Drop any sensitivity columns that were written in a previous run.
    mr_all <- mr_all[, setdiff(names(mr_all), setdiff(join_cols, c("cohort", "seqid"))),
                     with = FALSE]
    mr_all <- merge(mr_all, sens_all[, ..join_cols],
                    by = c("cohort", "seqid"), all.x = TRUE)
    fwrite(mr_all, mr_path, sep = "\t")
    n_pass <- sum(mr_all$passes_sensitivity, na.rm = TRUE)
    log_info("Enriched %s: %d/%d proteins pass sensitivity", mr_path, n_pass, nrow(mr_all))
  }
}

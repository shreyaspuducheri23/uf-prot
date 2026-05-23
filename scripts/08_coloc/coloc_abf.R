#!/usr/bin/env Rscript
# 08_coloc/coloc_abf.R
# Sensitivity colocalization using coloc.abf on ±1Mb cis regions.
#
# Usage:
#   Rscript scripts/08_coloc/coloc_abf.R [--cohort ARIC_EA] [--limit 50]

suppressPackageStartupMessages({
  library(coloc)
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
source(file.path(repo_root, "scripts", "rlib", "coloc_abf.R"))

setup_logger("08_coloc_abf")

# ── CLI args ──────────────────────────────────────────────────────────────────
args        <- commandArgs(trailingOnly = TRUE)
cohorts_all <- c("ARIC_EA", "deCODE", "UKB_PPP", "Fenland")
cohort_arg  <- "all"
limit_arg   <- Inf

i <- 1
while (i <= length(args)) {
  if (args[i] == "--cohort") { cohort_arg <- args[i+1]; i <- i + 2 }
  else if (args[i] == "--limit") { limit_arg <- as.integer(args[i+1]); i <- i + 2 }
  else i <- i + 1
}
run_cohorts <- if (cohort_arg == "all") cohorts_all else cohort_arg

KIM_N   <- 434152L
KIM_S   <- 74318 / (74318 + 359834)   # case fraction

run_cohort_coloc <- function(cohort) {
  region_base <- file.path("processed_data", "coloc", "regions", cohort)
  out_path    <- file.path("processed_data", "coloc", paste0("coloc_abf_", cohort, ".tsv"))
  state_path  <- file.path("processed_data", cohort, "_state_08_coloc_abf.rds")

  if (!dir.exists(region_base)) {
    log_warn("%s: no coloc regions directory", cohort)
    return(invisible(NULL))
  }

  candidates <- list.dirs(region_base, full.names = TRUE, recursive = FALSE)
  if (is.finite(limit_arg)) candidates <- head(candidates, limit_arg)

  cp <- checkpoint_load(state_path)
  todo <- candidates[basename(candidates) %in%
                      remaining(basename(candidates), cp)]

  log_info("%s: %d candidates, %d remaining", cohort, length(candidates), length(todo))

  results <- list()
  tick <- serial_progress(length(todo), label = paste(cohort, "coloc.abf"))

  for (region_dir in todo) {
    tick()
    seqid    <- basename(region_dir)
    exp_path <- file.path(region_dir, "exposure.tsv")
    out_fp   <- file.path(region_dir, "outcome.tsv")

    tryCatch({
      if (!file.exists(exp_path) || !file.exists(out_fp)) {
        cp <- checkpoint_mark(cp, seqid); next
      }

      exp_df <- fread(exp_path, sep = "\t",
                      colClasses = list(character = c("chrom", "rsid")))
      out_df <- fread(out_fp, sep = "\t",
                      colClasses = list(character = c("chromosome", "rsid")))

      # Build allele-aware keys — chrom:pos:EA:OA uniquely identifies each allele
      # combination. All exposure files are now in hg38 (lifted in steps 4/8a), so
      # positions are comparable with the hg38 outcome.
      exp_key <- paste(exp_df$chrom, exp_df$pos, exp_df$EA, exp_df$OA, sep = ":")

      out_key_fwd <- paste(out_df$chromosome, out_df$base_pair_location,
                           out_df$effect_allele, out_df$other_allele, sep = ":")
      out_key_rev <- paste(out_df$chromosome, out_df$base_pair_location,
                           out_df$other_allele, out_df$effect_allele, sep = ":")

      # Tag each outcome row with the exposure key it matches (forward or flipped).
      # fcase() short-circuits: forward match wins over reverse.
      out_df[, match_key := fcase(
        out_key_fwd %in% exp_key, out_key_fwd,
        out_key_rev %in% exp_key, out_key_rev,
        default = NA_character_
      )]
      # For reverse-matched rows, flip beta and EAF so effect direction is aligned
      out_df[match_key != out_key_fwd & !is.na(match_key), `:=`(
        beta                    = -beta,
        effect_allele_frequency =  1 - effect_allele_frequency
      )]
      out_df <- out_df[!is.na(match_key)]

      common <- intersect(
        exp_key[!duplicated(exp_key)],
        out_df$match_key[!duplicated(out_df$match_key)]
      )
      if (length(common) < 5) {
        cp <- checkpoint_mark(cp, seqid); next
      }

      exp_sub     <- exp_df[!duplicated(exp_key) & exp_key %in% common]
      exp_sub_key <- exp_key[!duplicated(exp_key) & exp_key %in% common]

      out_sub     <- out_df[!duplicated(out_df$match_key) & out_df$match_key %in% common]
      out_sub_key <- out_df$match_key[!duplicated(out_df$match_key) & out_df$match_key %in% common]

      # Align both frames by sorting on the shared key
      exp_sub <- exp_sub[order(exp_sub_key)]
      out_sub <- out_sub[order(out_sub_key)]
      exp_snp_key <- sort(exp_sub_key)

      # Infer N_exp from data
      N_exp <- if ("N" %in% names(exp_sub)) as.integer(median(exp_sub$N, na.rm = TRUE)) else 10000L

      coloc_res <- run_coloc_abf(
        exp_df = data.frame(
          snp    = exp_snp_key,
          beta   = as.numeric(exp_sub$beta),
          se     = as.numeric(exp_sub$se),
          eaf    = as.numeric(exp_sub$EAF)
        ),
        out_df = data.frame(
          snp    = sort(out_sub_key),
          beta   = as.numeric(out_sub$beta),
          se     = as.numeric(out_sub$standard_error),
          eaf    = as.numeric(out_sub$effect_allele_frequency)
        ),
        N_exp = N_exp,
        N_out = KIM_N
      )

      pp <- coloc_res$summary
      results[[seqid]] <- data.frame(
        cohort         = cohort,
        seqid          = seqid,
        n_snps         = length(common),
        PP_H0          = pp["PP.H0.abf"],
        PP_H1          = pp["PP.H1.abf"],
        PP_H2          = pp["PP.H2.abf"],
        PP_H3          = pp["PP.H3.abf"],
        PP_H4          = pp["PP.H4.abf"],
        coloc_positive = pp["PP.H4.abf"] >= 0.8,
        stringsAsFactors = FALSE
      )
      cp <- checkpoint_mark(cp, seqid)
    }, error = function(e) {
      log_warn("%s %s: coloc.abf failed — %s", cohort, seqid, conditionMessage(e))
      cp <<- checkpoint_mark(cp, seqid)
    })
  }

  if (length(results) == 0) return(invisible(NULL))

  dt <- rbindlist(results, fill = TRUE)
  if (file.exists(out_path)) {
    existing <- fread(out_path, sep = "\t")
    dt <- rbindlist(list(existing[!seqid %in% dt$seqid], dt), fill = TRUE)
  }
  fwrite(dt, out_path, sep = "\t")
  log_info("%s: %d coloc.abf results written (%d coloc_positive)",
           cohort, nrow(dt), sum(dt$coloc_positive, na.rm = TRUE))
  invisible(dt)
}

lapply(run_cohorts, run_cohort_coloc)

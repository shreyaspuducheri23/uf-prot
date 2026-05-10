# Wrappers around TwoSampleMR::harmonise_data.
# Called from Python via subprocess (Rscript) with TSV file paths as args.
#
# Usage (standalone):
#   Rscript scripts/rlib/harmonise.R --exp exposure.tsv --out outcome.tsv --result harmonised.tsv

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

harmonise_files <- function(exp_path, out_path, result_path) {
  exp_raw  <- fread(exp_path,  sep = "\t", colClasses = list(character = c("chrom", "rsid", "EA", "OA")))
  out_raw  <- fread(out_path,  sep = "\t", colClasses = list(character = c("chromosome", "rsid")))

  exp_fmt <- format_data(
    as.data.frame(exp_raw),
    type          = "exposure",
    snp_col       = "rsid",
    beta_col      = "beta",
    se_col        = "se",
    eaf_col       = "EAF",
    effect_allele_col = "EA",
    other_allele_col  = "OA",
    pval_col      = "pval",
    samplesize_col = "N",
    chr_col       = "chrom",
    pos_col       = "pos",
  )
  exp_fmt$exposure <- exp_raw$seqid[1]

  out_fmt <- format_data(
    as.data.frame(out_raw),
    type          = "outcome",
    snp_col       = "rsid",
    beta_col      = "beta",
    se_col        = "se",
    eaf_col       = "EAF",
    effect_allele_col = "EA",
    other_allele_col  = "OA",
    pval_col      = "pval",
    samplesize_col = "N",
    chr_col       = "chrom",
    pos_col       = "pos",
  )
  out_fmt$outcome <- "uterine_fibroids_kim2025"

  harmonised <- harmonise_data(exp_fmt, out_fmt, action = 2)
  fwrite(harmonised, result_path, sep = "\t")
  log_info("Harmonised %d SNPs → %s", nrow(harmonised), result_path)
}

if (!interactive()) {
  setup_logger("harmonise")
  args <- commandArgs(trailingOnly = TRUE)
  parser <- list(exp = NULL, out = NULL, result = NULL)
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--exp")    { parser$exp    <- args[i+1]; i <- i + 2 }
    else if (args[i] == "--out")    { parser$out    <- args[i+1]; i <- i + 2 }
    else if (args[i] == "--result") { parser$result <- args[i+1]; i <- i + 2 }
    else i <- i + 1
  }
  stopifnot(!is.null(parser$exp), !is.null(parser$out), !is.null(parser$result))
  harmonise_files(parser$exp, parser$out, parser$result)
}

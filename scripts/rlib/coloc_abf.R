# coloc.abf wrapper for sensitivity colocalization.
# Called from code/08_coloc/coloc_abf.R.

suppressPackageStartupMessages(library(coloc))

run_coloc_abf <- function(exp_df, out_df, N_exp, N_out, s_out, type_exp = "quant") {
  # exp_df / out_df: data frames with columns: snp, beta, se, eaf (EAF), N (optional)
  # Returns coloc.abf result list

  d1 <- list(
    beta   = exp_df$beta,
    varbeta = exp_df$se^2,
    snp    = exp_df$snp,
    MAF    = pmin(exp_df$eaf, 1 - exp_df$eaf),
    type   = type_exp,
    N      = N_exp
  )

  d2 <- list(
    beta    = out_df$beta,
    varbeta = out_df$se^2,
    snp     = out_df$snp,
    MAF     = pmin(out_df$eaf, 1 - out_df$eaf),
    type    = "cc",
    N       = N_out,
    s       = s_out
  )

  coloc.abf(d1, d2)
}

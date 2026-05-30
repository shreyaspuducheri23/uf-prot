# Shared MR method drivers: Wald, IVW-MRE, weighted median/mode, Egger, Steiger.
# Called from 06_mr/run_mr.R and 07_sensitivity/run_sensitivity.R.

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

run_protein_mr <- function(harm_dt) {
  # harm_dt: harmonised data.table for one protein (output of harmonise_data)
  n_snps <- sum(harm_dt$mr_keep)
  if (n_snps == 0) return(NULL)

  methods <- if (n_snps == 1) {
    c("mr_wald_ratio")
  } else if (n_snps == 2) {
    c("mr_ivw")
  } else {
    c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
  }

  res <- suppressMessages(mr(harm_dt, method_list = methods))
  if (nrow(res) == 0) return(NULL)

  # Primary estimate
  primary_method <- if (n_snps == 1) "Wald ratio" else "Inverse variance weighted"
  primary <- res[res$method == primary_method, ]
  if (nrow(primary) == 0) primary <- res[1, ]

  list(
    seqid       = harm_dt$exposure[1],
    n_snps      = n_snps,
    method      = primary$method,
    beta        = primary$b,
    se          = primary$se,
    pval        = primary$pval,
    OR          = exp(primary$b),
    OR_lo95     = exp(primary$b - 1.96 * primary$se),
    OR_hi95     = exp(primary$b + 1.96 * primary$se),
    all_results = res
  )
}


run_sensitivity <- function(harm_dt) {
  n_snps <- sum(harm_dt$mr_keep)
  out <- list(seqid = harm_dt$exposure[1], n_snps = n_snps)

  if (n_snps >= 2) {
    het <- mr_heterogeneity(harm_dt)
    ivw_het <- het[het$method == "Inverse variance weighted", ]
    out$Q          <- if (nrow(ivw_het) > 0) ivw_het$Q       else NA
    out$Q_df       <- if (nrow(ivw_het) > 0) ivw_het$Q_df    else NA
    out$Q_pval     <- if (nrow(ivw_het) > 0) ivw_het$Q_pval  else NA
    out$I2         <- if (!is.null(out$Q) && !is.null(out$Q_df) && out$Q_df > 0)
                        max(0, (out$Q - out$Q_df) / out$Q) else NA
    out$heterogeneous <- !is.na(out$I2) && !is.na(out$Q_pval) &&
                          out$I2 >= 0.5 && out$Q_pval < 0.05
  }

  if (n_snps >= 3) {
    plei <- mr_pleiotropy_test(harm_dt)
    out$egger_intercept       <- if (nrow(plei) > 0) plei$egger_intercept  else NA
    out$egger_intercept_se    <- if (nrow(plei) > 0) plei$se               else NA
    out$egger_intercept_pval  <- if (nrow(plei) > 0) plei$pval             else NA
    out$directional_pleiotropy <- !is.na(out$egger_intercept_pval) &&
                                   out$egger_intercept_pval < 0.05

    # Alternative methods
    res_all <- suppressMessages(mr(harm_dt, method_list = c("mr_ivw", "mr_weighted_median",
                                                             "mr_weighted_mode", "mr_egger_regression")))
    # Direction consistency across all methods
    signs <- sign(res_all$b)
    out$direction_consistent <- length(unique(signs)) == 1
  }

  # Steiger filtering
  steiger <- tryCatch(directionality_test(harm_dt), error = function(e) NULL)
  if (!is.null(steiger)) {
    out$steiger_correct  <- steiger$correct_causal_direction
    out$steiger_pval     <- steiger$steiger_pval
  }

  out$passes_sensitivity <- if (n_snps == 1) {
    isTRUE(out$steiger_correct != FALSE)
  } else if (n_snps == 2) {
    isTRUE(out$steiger_correct != FALSE) &&
      !isTRUE(out$heterogeneous)
  } else {
    isTRUE(out$steiger_correct != FALSE) &&
      !isTRUE(out$heterogeneous) &&
      !isTRUE(out$directional_pleiotropy) &&
      !identical(out$direction_consistent, FALSE)
  }

  out
}


bh_fdr <- function(pvals, alpha = 0.05) {
  p.adjust(pvals, method = "BH")
}

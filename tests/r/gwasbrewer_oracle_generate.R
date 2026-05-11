#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GWASBrewer)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript gwasbrewer_oracle_generate.R <out_dir> <scenario> <seed> <n_proteins> <snps_per_protein> [locus_spacing_bp]")
}

out_dir <- args[[1]]
scenario <- args[[2]]
seed <- as.integer(args[[3]])
n_proteins <- as.integer(args[[4]])
snps_per_protein <- as.integer(args[[5]])
locus_spacing_bp <- if (length(args) >= 6) as.integer(args[[6]]) else 2000000L

if (is.na(seed) || is.na(n_proteins) || is.na(snps_per_protein) || is.na(locus_spacing_bp)) {
  stop("seed, n_proteins, snps_per_protein, and locus_spacing_bp must be integers")
}

if (n_proteins < 8) {
  stop("n_proteins must be >= 8 so all oracle classes are represented")
}

set.seed(seed)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

classes <- rep("null", n_proteins)
classes[1:4] <- "causal"
classes[5:6] <- "pleiotropic"
classes[7:8] <- "weak"

mediation_scale <- switch(
  scenario,
  shared_signal = 0.85,
  coloc_mismatch = 0.85,
  pleiotropy_guardrail = 0.85,
  weak_instrument_guardrail = 0.70,
  proxy_branch = 0.90,
  full_pipeline = 0.80,
  null_control = 0.0,
  0.85
)
pleiotropy_beta <- switch(
  scenario,
  pleiotropy_guardrail = 0.28,
  proxy_branch = 0.18,
  full_pipeline = 0.15,
  0.20
)

exposure_rows <- list()
outcome_rows <- list()
manifest_rows <- list()

for (i in seq_len(n_proteins)) {
  seqid <- sprintf("SeqId_%02d", i - 1)
  gene <- sprintf("G%02d", i - 1)
  uniprot <- sprintf("P%05d", i - 1)
  class_i <- classes[[i]]
  chrom <- "1"
  tss <- 1000000L + (i - 1L) * locus_spacing_bp

  J <- snps_per_protein
  G <- matrix(c(0, mediation_scale, 0, 0), nrow = 2, byrow = TRUE)
  if (class_i == "null") {
    G <- matrix(0, nrow = 2, ncol = 2)
  }

  sim <- sim_mv(
    N = c(7000, 434152),
    J = J,
    h2 = c(0.20, 0.20),
    pi = c(0.30, 0.30),
    G = G,
    sporadic_pleiotropy = FALSE,
    pi_exact = TRUE
  )

  beta_exp <- as.numeric(sim$beta_hat[, 1])
  se_exp <- as.numeric(sim$se_beta_hat[, 1])

  if (class_i == "weak") {
    beta_exp <- beta_exp * 0.12
    se_exp <- pmax(se_exp * 2.0, 0.04)
  }

  z_exp <- abs(beta_exp) / pmax(se_exp, 1e-12)
  pval_exp <- pmax(2 * pnorm(-z_exp), 1e-300)
  if (class_i %in% c("weak", "null")) {
    pval_exp <- pmax(pval_exp, 1e-7)
  } else {
    n_sig <- min(4, J)
    pval_exp[seq_len(n_sig)] <- 1e-10
    if (J > n_sig) {
      pval_exp[(n_sig + 1):J] <- pmax(pval_exp[(n_sig + 1):J], 1e-7)
    }
  }

  beta_out <- as.numeric(sim$beta_hat[, 2])
  se_out <- as.numeric(sim$se_beta_hat[, 2])

  if (class_i %in% c("causal", "weak")) {
    beta_out <- mediation_scale * beta_exp + rnorm(J, mean = 0, sd = 0.02)
  }
  if (class_i == "pleiotropic") {
    beta_out <- pleiotropy_beta + rnorm(J, mean = 0, sd = 0.02)
  }
  if (class_i == "null") {
    beta_out <- rnorm(J, mean = 0, sd = 0.02)
  }

  z_out <- abs(beta_out) / pmax(se_out, 1e-12)
  pval_out <- pmax(2 * pnorm(-z_out), 1e-300)

  ea_vec <- rep("A", J)
  oa_vec <- rep("G", J)
  eaf_vec <- pmin(pmax(runif(J, 0.08, 0.42), 0.01), 0.99)

  mean_f <- mean((beta_exp / pmax(se_exp, 1e-12))^2)
  w <- 1 / pmax(se_out, 1e-12)^2
  denom <- sum(w * beta_exp^2)
  ivw_beta <- ifelse(denom > 0, sum(w * beta_exp * beta_out) / denom, 0.0)
  target_beta <- if (class_i == "causal") ivw_beta else 0.0

  expected_sign <- if (target_beta > 0) 1L else if (target_beta < 0) -1L else 0L
  expected_tier <- if (class_i == "causal" && scenario != "coloc_mismatch") {
    "Tier1_replicated"
  } else if (class_i == "causal" && scenario == "coloc_mismatch") {
    "Tier2"
  } else {
    "Tier3"
  }

  manifest_rows[[i]] <- data.table(
    seqid = seqid,
    gene = gene,
    uniprot = uniprot,
    class = class_i,
    chrom = chrom,
    tss = tss,
    expected_effect_sign = expected_sign,
    target_mr_beta = target_beta,
    beta_tolerance = ifelse(class_i == "causal", 0.08, 0.12),
    expected_tier_tendency = expected_tier,
    expected_coloc_positive = as.logical(class_i == "causal" && scenario != "coloc_mismatch"),
    expected_weak = as.logical(class_i == "weak"),
    expected_null = as.logical(class_i == "null"),
    expected_proxy = as.logical(scenario == "proxy_branch" && i %% 2 == 0),
    mediation_scale = mediation_scale,
    pleiotropy_beta = pleiotropy_beta,
    mean_F_expected = mean_f
  )

  for (j in seq_len(J)) {
    pos <- tss - 350000L + (j - 1L) * 100000L
    rsid <- sprintf("rs%s_%02d", seqid, j - 1)

    exposure_rows[[length(exposure_rows) + 1L]] <- data.table(
      seqid = seqid,
      gene = gene,
      uniprot = uniprot,
      class = class_i,
      chrom = chrom,
      pos = pos,
      rsid = rsid,
      EA = ea_vec[[j]],
      OA = oa_vec[[j]],
      EAF = eaf_vec[[j]],
      beta = beta_exp[[j]],
      se = se_exp[[j]],
      pval = pval_exp[[j]],
      N = 7000L,
      build = "hg19"
    )

    outcome_rows[[length(outcome_rows) + 1L]] <- data.table(
      seqid = seqid,
      class = class_i,
      chromosome = chrom,
      base_pair_location = pos,
      effect_allele = ea_vec[[j]],
      other_allele = oa_vec[[j]],
      beta = beta_out[[j]],
      standard_error = se_out[[j]],
      effect_allele_frequency = eaf_vec[[j]],
      p_value = pval_out[[j]],
      rsid = rsid,
      rs_id = rsid,
      hm_coordinate_conversion = "",
      hm_code = "",
      variant_id = ""
    )
  }
}

exp_dt <- rbindlist(exposure_rows, use.names = TRUE, fill = TRUE)
out_dt <- rbindlist(outcome_rows, use.names = TRUE, fill = TRUE)
manifest_dt <- rbindlist(manifest_rows, use.names = TRUE, fill = TRUE)

fwrite(exp_dt, file.path(out_dir, "exposure_gwas.tsv"), sep = "\t")
fwrite(out_dt, file.path(out_dir, "outcome_gwas.tsv"), sep = "\t")
fwrite(manifest_dt, file.path(out_dir, "oracle_manifest.tsv"), sep = "\t")

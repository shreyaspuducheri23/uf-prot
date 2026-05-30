# Allele-aware summary-statistic alignment for coloc inputs.

align_coloc_region <- function(exp_df, out_df) {
  exp_df <- as.data.frame(exp_df)
  out_df <- as.data.frame(out_df)

  exp_key <- paste(exp_df$chrom, exp_df$pos, exp_df$EA, exp_df$OA, sep = ":")

  out_key_fwd <- paste(out_df$chromosome, out_df$base_pair_location,
                       out_df$effect_allele, out_df$other_allele, sep = ":")
  out_key_rev <- paste(out_df$chromosome, out_df$base_pair_location,
                       out_df$other_allele, out_df$effect_allele, sep = ":")

  exp_keys <- unique(exp_key)
  match_key <- rep(NA_character_, length(out_key_fwd))
  fwd <- out_key_fwd %in% exp_keys
  rev <- !fwd & out_key_rev %in% exp_keys
  match_key[fwd] <- out_key_fwd[fwd]
  match_key[rev] <- out_key_rev[rev]

  out_aligned_all <- out_df[!is.na(match_key), , drop = FALSE]
  out_match_key <- match_key[!is.na(match_key)]
  if (nrow(out_aligned_all) == 0) {
    return(list(
      exp = exp_df[0, , drop = FALSE],
      out = out_df[0, , drop = FALSE],
      snp = character()
    ))
  }

  rev_kept <- rev[!is.na(match_key)]
  if (any(rev_kept)) {
    out_aligned_all$beta[rev_kept] <- -as.numeric(out_aligned_all$beta[rev_kept])
    out_aligned_all$effect_allele_frequency[rev_kept] <-
      1 - as.numeric(out_aligned_all$effect_allele_frequency[rev_kept])
  }
  out_aligned_all$match_key <- out_match_key

  exp_first <- !duplicated(exp_key)
  out_first <- !duplicated(out_aligned_all$match_key)
  common <- intersect(exp_key[exp_first], out_aligned_all$match_key[out_first])
  if (length(common) == 0) {
    return(list(
      exp = exp_df[0, , drop = FALSE],
      out = out_df[0, , drop = FALSE],
      snp = character()
    ))
  }

  exp_sub <- exp_df[exp_first & exp_key %in% common, , drop = FALSE]
  exp_sub_key <- exp_key[exp_first & exp_key %in% common]
  out_sub <- out_aligned_all[out_first & out_aligned_all$match_key %in% common, , drop = FALSE]
  out_sub_key <- out_aligned_all$match_key[out_first & out_aligned_all$match_key %in% common]

  exp_sub <- exp_sub[order(exp_sub_key), , drop = FALSE]
  out_sub <- out_sub[order(out_sub_key), , drop = FALSE]
  out_sub$match_key <- NULL

  list(exp = exp_sub, out = out_sub, snp = sort(common))
}

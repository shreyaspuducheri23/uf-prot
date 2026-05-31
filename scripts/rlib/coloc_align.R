# Allele-aware summary-statistic alignment for coloc inputs.

normalize_allele <- function(allele) {
  if (is.logical(allele)) {
    normalized <- rep(NA_character_, length(allele))
    normalized[!is.na(allele) & allele] <- "T"
    normalized[!is.na(allele) & !allele] <- "F"
    return(normalized)
  }

  toupper(as.character(allele))
}

complement_allele <- function(allele) {
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  unname(comp[normalize_allele(allele)])
}

align_coloc_region <- function(exp_df, out_df) {
  exp_df <- as.data.frame(exp_df)
  out_df <- as.data.frame(out_df)

  exp_effect <- normalize_allele(exp_df$EA)
  exp_other <- normalize_allele(exp_df$OA)
  out_effect <- normalize_allele(out_df$effect_allele)
  out_other <- normalize_allele(out_df$other_allele)

  exp_key <- paste(exp_df$chrom, exp_df$pos, exp_effect, exp_other, sep = ":")

  out_key_fwd <- paste(out_df$chromosome, out_df$base_pair_location,
                       out_effect, out_other, sep = ":")
  out_key_rev <- paste(out_df$chromosome, out_df$base_pair_location,
                       out_other, out_effect, sep = ":")
  out_effect_comp <- complement_allele(out_effect)
  out_other_comp <- complement_allele(out_other)
  out_has_complement <- !is.na(out_effect_comp) & !is.na(out_other_comp)
  out_key_comp_fwd <- paste(out_df$chromosome, out_df$base_pair_location,
                            out_effect_comp, out_other_comp, sep = ":")
  out_key_comp_rev <- paste(out_df$chromosome, out_df$base_pair_location,
                            out_other_comp, out_effect_comp, sep = ":")

  exp_keys <- unique(exp_key)
  match_key <- rep(NA_character_, length(out_key_fwd))
  fwd <- out_key_fwd %in% exp_keys
  rev <- !fwd & out_key_rev %in% exp_keys
  comp_fwd <- !fwd & !rev & out_has_complement & out_key_comp_fwd %in% exp_keys
  comp_rev <- !fwd & !rev & !comp_fwd & out_has_complement & out_key_comp_rev %in% exp_keys
  match_key[fwd] <- out_key_fwd[fwd]
  match_key[rev] <- out_key_rev[rev]
  match_key[comp_fwd] <- out_key_comp_fwd[comp_fwd]
  match_key[comp_rev] <- out_key_comp_rev[comp_rev]

  out_aligned_all <- out_df[!is.na(match_key), , drop = FALSE]
  out_match_key <- match_key[!is.na(match_key)]
  if (nrow(out_aligned_all) == 0) {
    return(list(
      exp = exp_df[0, , drop = FALSE],
      out = out_df[0, , drop = FALSE],
      snp = character()
    ))
  }

  flip_kept <- (rev | comp_rev)[!is.na(match_key)]
  if (any(flip_kept)) {
    out_aligned_all$beta[flip_kept] <- -as.numeric(out_aligned_all$beta[flip_kept])
    out_aligned_all$effect_allele_frequency[flip_kept] <-
      1 - as.numeric(out_aligned_all$effect_allele_frequency[flip_kept])
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

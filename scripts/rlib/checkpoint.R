# Per-unit RDS-based checkpointing for R scripts.
# Usage:
#   cp <- checkpoint_load("processed_data/ARIC_EA/_state_06.rds")
#   for (seqid in remaining(all_seqids, cp)) {
#     ... process ...
#     checkpoint_mark(cp, seqid)
#   }

checkpoint_load <- function(state_path) {
  if (file.exists(state_path)) {
    tryCatch(readRDS(state_path), error = function(e) list(done = character(0), path = state_path))
  } else {
    list(done = character(0), path = state_path)
  }
}

checkpoint_mark <- function(cp, key) {
  cp$done <- unique(c(cp$done, key))
  dir.create(dirname(cp$path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(cp, cp$path)
  invisible(cp)
}

remaining <- function(keys, cp) {
  keys[!keys %in% cp$done]
}

checkpoint_n_done <- function(cp) length(cp$done)

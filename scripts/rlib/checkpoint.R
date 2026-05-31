# Per-unit RDS-based checkpointing for R scripts.
# Usage:
#   cp <- checkpoint_load("processed_data/ARIC_EA/_state_06.rds")
#   for (seqid in remaining(all_seqids, cp)) {
#     ... process ...
#     checkpoint_mark(cp, seqid)
#   }

checkpoint_load <- function(state_path) {
  if (file.exists(state_path)) {
    cp <- tryCatch(readRDS(state_path), error = function(e) list(done = character(0), path = state_path))
  } else {
    cp <- list(done = character(0), path = state_path)
  }
  if (is.null(cp$done)) cp$done <- character(0)
  if (is.null(cp$failed)) cp$failed <- list()
  cp$path <- state_path
  cp
}

checkpoint_mark <- function(cp, key) {
  cp$done <- unique(c(cp$done, key))
  if (!is.null(cp$failed[[key]])) cp$failed[[key]] <- NULL
  dir.create(dirname(cp$path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(cp, cp$path)
  invisible(cp)
}

checkpoint_mark_failed <- function(cp, key, reason = "") {
  cp$done <- setdiff(cp$done, key)
  cp$failed[[key]] <- list(
    reason = as.character(reason),
    updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z")
  )
  dir.create(dirname(cp$path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(cp, cp$path)
  invisible(cp)
}

remaining <- function(keys, cp, include_failed = TRUE) {
  if (include_failed) {
    keys[!keys %in% cp$done]
  } else {
    keys[!keys %in% cp$done & !keys %in% names(cp$failed)]
  }
}

checkpoint_n_done <- function(cp) length(cp$done)
checkpoint_n_failed <- function(cp) length(cp$failed)

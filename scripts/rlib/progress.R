# Progress tracking for R scripts using progressr + pbmcapply.

# For parallel lapply (returns list in same order):
# results <- pbmc_lapply(items, function(x) { ... }, n_workers = 4)
#
# For serial progress:
# p <- serial_progress(n_items, "Processing proteins")
# for (...) { p(); ... }

pbmc_lapply <- function(items, fn, n_workers = 1L, ...) {
  if (!requireNamespace("pbmcapply", quietly = TRUE)) {
    stop("Install pbmcapply: install.packages('pbmcapply')")
  }
  pbmcapply::pbmclapply(items, fn, mc.cores = n_workers, ...)
}

serial_progress <- function(n, label = "") {
  if (requireNamespace("progress", quietly = TRUE)) {
    pb <- progress::progress_bar$new(
      format = paste0(label, " [:bar] :current/:total (:percent) ETA :eta"),
      total = n, clear = FALSE
    )
    function() pb$tick()
  } else {
    # Fallback: print every 10%
    i <- 0L
    step <- max(1L, as.integer(n / 10))
    function() {
      i <<- i + 1L
      if (i %% step == 0L) cat(sprintf("%s: %d/%d\n", label, i, n))
    }
  }
}

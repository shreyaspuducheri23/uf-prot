# Shared logging setup for R scripts: ISO-8601 timestamps + dual sink to logs/
# Usage: source("scripts/rlib/logging.R"); setup_logger("06_mr")

library(futile.logger)

setup_logger <- function(step_name, log_dir = "logs") {
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
  log_path <- file.path(log_dir, paste0(step_name, "_", ts, ".log"))

  flog.appender(appender.tee(log_path), name = "ROOT")
  flog.layout(layout.format(paste0("[~t] ~l ~m")), name = "ROOT")
  flog.threshold(INFO, name = "ROOT")

  flog.info("Log file: %s", log_path)
  invisible(log_path)
}

log_info  <- function(...) flog.info(...,  name = "ROOT")
log_warn  <- function(...) flog.warn(...,  name = "ROOT")
log_error <- function(...) flog.error(..., name = "ROOT")

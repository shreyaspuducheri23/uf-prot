# Shared pipeline configuration for R entry points.

.pipeline_config_env <- new.env(parent = emptyenv())

pipeline_config_path <- function() {
  root <- if (exists("repo_root", inherits = TRUE)) {
    get("repo_root", inherits = TRUE)
  } else {
    normalizePath(".")
  }
  file.path(root, "config", "pipeline.json")
}

pipeline_config <- function(path = pipeline_config_path()) {
  path <- normalizePath(path, mustWork = TRUE)
  if (!exists("path", envir = .pipeline_config_env, inherits = FALSE) ||
      !identical(.pipeline_config_env$path, path)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("R package 'jsonlite' is required to read config/pipeline.json", call. = FALSE)
    }
    .pipeline_config_env$path <- path
    .pipeline_config_env$config <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  }
  .pipeline_config_env$config
}

pipeline_cohorts <- function() {
  names(pipeline_config()$cohorts)
}

pipeline_kim_n <- function() {
  as.integer(pipeline_config()$outcome$kim_N)
}

pipeline_kim_s <- function() {
  outcome <- pipeline_config()$outcome
  cases <- as.numeric(outcome$kim_cases)
  controls <- as.numeric(outcome$kim_controls)
  cases / (cases + controls)
}

pipeline_kim_ncase <- function() {
  as.numeric(pipeline_config()$outcome$kim_cases)
}

pipeline_kim_ncontrol <- function() {
  as.numeric(pipeline_config()$outcome$kim_controls)
}

pipeline_prevalence <- function() {
  as.numeric(pipeline_config()$outcome$prevalence)
}

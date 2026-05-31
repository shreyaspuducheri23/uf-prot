"""Regression tests for scripts/rlib/mr_methods.R method selection."""
import json
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def _require_r_twosamplemr() -> None:
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not available")
    probe = subprocess.run(
        [
            "Rscript",
            "-e",
            'quit(save="no", status=ifelse(requireNamespace("TwoSampleMR", quietly=TRUE),0,42))',
        ],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    if probe.returncode == 42:
        pytest.skip("TwoSampleMR is not installed")
    assert probe.returncode == 0, probe.stderr


def test_mr_methods_use_explicit_ivw_mre(tmp_path):
    _require_r_twosamplemr()

    script = tmp_path / "check_mr_methods.R"
    script.write_text(
        textwrap.dedent(
            f"""
            repo_root <- {json.dumps(str(REPO_ROOT))}
            source(file.path(repo_root, "scripts", "rlib", "mr_methods.R"))

            calls <- list()
            make_harm <- function(n) {{
              data.table::data.table(
                exposure = rep(paste0("SeqId_", n), n),
                mr_keep = rep(TRUE, n)
              )
            }}

            mr <- function(dat, method_list) {{
              calls[[length(calls) + 1]] <<- method_list
              if (identical(method_list, c("mr_wald_ratio"))) {{
                return(data.frame(
                  method = "Wald ratio",
                  b = 0.2,
                  se = 0.1,
                  pval = 0.04,
                  stringsAsFactors = FALSE
                ))
              }}
              data.frame(
                method = c(
                  "Inverse variance weighted",
                  "Inverse variance weighted (multiplicative random effects)",
                  "Weighted median",
                  "Weighted mode",
                  "MR Egger"
                ),
                b = c(9, 0.3, 0.3, 0.3, 0.3),
                se = c(9, 0.2, 0.2, 0.2, 0.2),
                pval = c(1, 0.05, 0.05, 0.05, 0.05),
                stringsAsFactors = FALSE
              )
            }}

            res1 <- run_protein_mr(make_harm(1))
            stopifnot(identical(calls[[1]], c("mr_wald_ratio")))
            stopifnot(identical(as.character(res1$method), "Wald ratio"))

            res2 <- run_protein_mr(make_harm(2))
            stopifnot(identical(calls[[2]], c("mr_ivw_mre")))
            stopifnot(identical(
              as.character(res2$method),
              "Inverse variance weighted (multiplicative random effects)"
            ))
            stopifnot(identical(as.numeric(res2$beta), 0.3))

            res3 <- run_protein_mr(make_harm(3))
            stopifnot(identical(
              calls[[3]],
              c("mr_ivw_mre", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
            ))
            stopifnot(identical(
              as.character(res3$method),
              "Inverse variance weighted (multiplicative random effects)"
            ))

            calls <- list()
            mr_heterogeneity <- function(dat) {{
              data.frame(
                method = "Inverse variance weighted",
                Q = 1,
                Q_df = 2,
                Q_pval = 0.6,
                stringsAsFactors = FALSE
              )
            }}
            mr_pleiotropy_test <- function(dat) {{
              data.frame(egger_intercept = 0, se = 1, pval = 0.9)
            }}
            directionality_test <- function(dat) {{
              data.frame(correct_causal_direction = TRUE, steiger_pval = 0.5)
            }}

            sens <- run_sensitivity(make_harm(3))
            stopifnot(identical(
              calls[[1]],
              c("mr_ivw_mre", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
            ))
            stopifnot(isTRUE(sens$direction_consistent))
            stopifnot(isTRUE(sens$passes_sensitivity))
            """
        )
    )

    result = subprocess.run(
        ["Rscript", str(script)],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert result.returncode == 0, f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

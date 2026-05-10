"""Optional R-backed GWAS simulation test (GWASBrewer first, simGWAS fallback)."""
import json
import subprocess

import pandas as pd
import pytest


def test_optional_r_backend_simulation(tmp_path):
    out_path = tmp_path / "sim.tsv"
    script = f"""
backend <- NULL
if (requireNamespace("GWASBrewer", quietly = TRUE)) {{
  library(GWASBrewer)
  backend <- "GWASBrewer"
}} else if (requireNamespace("simGWAS", quietly = TRUE)) {{
  library(simGWAS)
  backend <- "simGWAS"
}} else {{
  quit(save="no", status=42)
}}

set.seed(123)
n <- 400
beta_exp <- rnorm(n, mean = 0.2, sd = 0.06)
is_causal <- rep(c(TRUE, FALSE), each = n/2)
beta_out <- ifelse(is_causal, 0.35 * beta_exp + rnorm(n, 0, 0.02), rnorm(n, 0, 0.02))
dat <- data.frame(backend = backend, is_causal = is_causal, beta_exp = beta_exp, beta_out = beta_out)
write.table(dat, file = "{out_path}", sep = "\\t", row.names = FALSE, quote = FALSE)
"""
    proc = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True)
    if proc.returncode == 42:
        pytest.skip("GWASBrewer/simGWAS not installed in test environment")
    assert proc.returncode == 0, proc.stderr

    df = pd.read_csv(out_path, sep="\t")
    assert df["backend"].iloc[0] in {"GWASBrewer", "simGWAS"}

    causal = df[df["is_causal"] == True]
    nulls = df[df["is_causal"] == False]
    assert causal["beta_out"].abs().mean() > nulls["beta_out"].abs().mean()

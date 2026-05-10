"""Subprocess wrappers for PLINK 1.9 and PLINK2."""
import logging
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from scripts.lib.paths import LD_REF_PREFIX, PLINK1

log = logging.getLogger(__name__)

_PLINK1 = str(PLINK1)


def _run(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess:
    log.debug("plink cmd: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    if result.returncode != 0:
        log.error("plink stderr:\n%s", result.stderr)
        raise RuntimeError(f"plink exited with code {result.returncode}")
    return result


def clump(sumstats: pd.DataFrame, seqid: str,
          bfile: Path = LD_REF_PREFIX,
          window_kb: int = 1000, r2: float = 0.001, p1: float = 5e-8,
          pval_col: str = "pval", snp_col: str = "rsid",
          ) -> pd.DataFrame:
    """
    LD clumping via plink1.9.
    Returns a DataFrame of independent lead SNPs (subset of input rows).
    Drops rows where rsid == '.'.
    """
    df = sumstats[sumstats[snp_col] != "."].copy()
    if df.empty:
        return df

    with tempfile.TemporaryDirectory(prefix=f"clump_{seqid}_") as tmpdir:
        tmp = Path(tmpdir)
        # plink --clump needs SNP + P columns
        assoc = tmp / "assoc.tsv"
        df[[snp_col, pval_col]].rename(columns={snp_col: "SNP", pval_col: "P"}).to_csv(
            assoc, sep="\t", index=False
        )

        out_prefix = tmp / "clumped"
        cmd = [
            _PLINK1,
            "--bfile", str(bfile),
            "--clump", str(assoc),
            "--clump-kb", str(window_kb),
            "--clump-r2", str(r2),
            "--clump-p1", str(p1),
            "--clump-snp-field", "SNP",
            "--clump-field", "P",
            "--out", str(out_prefix),
        ]
        _run(cmd)

        clump_file = out_prefix.with_suffix(".clumped")
        if not clump_file.exists():
            log.warning(f"No clumped output for {seqid} (no variants passed threshold)")
            return pd.DataFrame(columns=df.columns)

        clumped = pd.read_csv(clump_file, sep=r"\s+")
        lead_rsids = set(clumped["SNP"].astype(str))

    return df[df[snp_col].isin(lead_rsids)].copy()


def find_proxies(missing_rsids: list[str],
                 bfile: Path = LD_REF_PREFIX,
                 ld_window_kb: int = 5000, r2_threshold: float = 0.8,
                 ) -> dict[str, tuple[str, float]]:
    """
    For each missing rsID, find the best proxy in bfile with r² >= r2_threshold.
    Returns {target_rsid: (proxy_rsid, r2)} for those that have a proxy.
    """
    if not missing_rsids:
        return {}

    with tempfile.TemporaryDirectory(prefix="proxies_") as tmpdir:
        tmp = Path(tmpdir)
        snplist = tmp / "targets.txt"
        snplist.write_text("\n".join(missing_rsids))

        out_prefix = tmp / "ld"
        cmd = [
            _PLINK1,
            "--bfile", str(bfile),
            "--show-tags", str(snplist),
            "--tag-r2", str(r2_threshold),
            "--tag-kb", str(ld_window_kb),
            "--list-all",
            "--out", str(out_prefix),
        ]
        _run(cmd)

        tags_file = out_prefix.with_suffix(".tags.list")
        if not tags_file.exists():
            return {}

        proxies: dict[str, tuple[str, float]] = {}
        with open(tags_file) as fh:
            header = fh.readline().split()
            for line in fh:
                parts = line.split()
                if len(parts) < 4:
                    continue
                target = parts[0]
                tags_str = parts[3]  # "proxy1|r2,proxy2|r2,..."
                best_rsid, best_r2 = None, 0.0
                for tag_entry in tags_str.split(","):
                    if "|" not in tag_entry:
                        continue
                    tag_rsid, tag_r2_str = tag_entry.rsplit("|", 1)
                    try:
                        tag_r2 = float(tag_r2_str)
                    except ValueError:
                        continue
                    if tag_rsid != target and tag_r2 > best_r2:
                        best_r2 = tag_r2
                        best_rsid = tag_rsid
                if best_rsid:
                    proxies[target] = (best_rsid, best_r2)

    return proxies


def r_square_matrix(snp_list: list[str], bfile: Path = LD_REF_PREFIX,
                    ) -> pd.DataFrame:
    """
    Compute pairwise r² matrix for snp_list via plink2 --r square.
    Returns a DataFrame (index and columns = rsIDs).
    """
    with tempfile.TemporaryDirectory(prefix="ld_matrix_") as tmpdir:
        tmp = Path(tmpdir)
        snpfile = tmp / "snps.txt"
        snpfile.write_text("\n".join(snp_list))

        out_prefix = tmp / "ld"
        # plink2 required for --r square; use system plink2 if available
        cmd = [
            "plink2",
            "--bfile", str(bfile),
            "--extract", str(snpfile),
            "--r", "square",
            "--out", str(out_prefix),
        ]
        _run(cmd)

        matrix_file = out_prefix.with_suffix(".r.square")
        if not matrix_file.exists():
            raise FileNotFoundError(f"plink2 r-square output not found: {matrix_file}")

        mat = pd.read_csv(matrix_file, sep="\t", header=None)
        mat.index = snp_list
        mat.columns = snp_list
        return mat

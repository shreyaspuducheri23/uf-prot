"""Subprocess wrappers for PLINK2."""
import logging
import re
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path

import pandas as pd

from scripts.lib.paths import LD_REF_PREFIX, PLINK2

log = logging.getLogger(__name__)

_PLINK2 = str(PLINK2)


@lru_cache(maxsize=4)
def _bim_pos_to_rsid(bfile: Path) -> dict[tuple[str, int], str]:
    """Load bfile .bim and return {(chrom, pos): rsid} for real rsIDs only."""
    bim = pd.read_csv(
        Path(str(bfile) + ".bim"),
        sep="\t",
        header=None,
        names=["chrom", "rsid", "cm", "pos", "a1", "a2"],
        dtype={"chrom": str, "rsid": str, "pos": int},
        usecols=["chrom", "rsid", "pos"],
    )
    bim = bim[bim["rsid"].str.startswith("rs")]
    return {(row.chrom, row.pos): row.rsid for row in bim.itertuples(index=False)}


def _run(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess:
    log.debug("plink cmd: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    if result.returncode != 0:
        log.error("plink stderr:\n%s", result.stderr)
        raise RuntimeError(f"plink exited with code {result.returncode}")
    return result


def _pick_r2_col(columns: list[str]) -> str | None:
    candidates = ["UNPHASED_R2", "PHASED_R2", "R2"]
    for col in candidates:
        if col in columns:
            return col
    return next((col for col in columns if col.endswith("_R2")), None)


def clump(
    sumstats: pd.DataFrame,
    seqid: str,
    bfile: Path = LD_REF_PREFIX,
    window_kb: int = 1000,
    r2: float = 0.001,
    p1: float = 5e-8,
    pval_col: str = "pval",
    snp_col: str = "rsid",
) -> pd.DataFrame:
    """
    LD clumping via plink2.
    Returns a DataFrame of independent lead SNPs (subset of input rows).
    For variants with rsid == '.', attempts position-based rsID annotation from
    the bfile .bim before clumping. Variants that cannot be matched are dropped.
    """
    df = sumstats.copy()
    missing_mask = df[snp_col] == "."
    if missing_mask.any() and "chrom" in df.columns and "pos" in df.columns:
        pos_map = _bim_pos_to_rsid(bfile)
        df.loc[missing_mask, snp_col] = df.loc[missing_mask].apply(
            lambda r: pos_map.get((str(r["chrom"]), int(r["pos"])), "."), axis=1
        )
        n_annotated = (df[snp_col] != ".").sum() - (~missing_mask).sum()
        log.debug(f"{seqid}: annotated {n_annotated}/{missing_mask.sum()} missing rsIDs from bim")
    df = df[df[snp_col] != "."].copy()
    if df.empty:
        return df

    with tempfile.TemporaryDirectory(prefix=f"clump_{seqid}_") as tmpdir:
        tmp = Path(tmpdir)
        assoc = tmp / "assoc.tsv"
        df[[snp_col, pval_col]].rename(columns={snp_col: "SNP", pval_col: "P"}).to_csv(
            assoc, sep="\t", index=False
        )

        out_prefix = tmp / "clumped"
        cmd = [
            _PLINK2,
            "--bfile",
            str(bfile),
            "--clump",
            str(assoc),
            "--clump-kb",
            str(window_kb),
            "--clump-r2",
            str(r2),
            "--clump-p1",
            str(p1),
            "--clump-id-field",
            "SNP",
            "--clump-p-field",
            "P",
            "--out",
            str(out_prefix),
        ]
        _run(cmd)

        clump_file = out_prefix.with_suffix(".clumps")
        if not clump_file.exists():
            log.warning(f"No clump output for {seqid} (no variants passed threshold)")
            return pd.DataFrame(columns=df.columns)

        clumped = pd.read_csv(clump_file, sep=r"\s+")
        snp_col_name = "ID" if "ID" in clumped.columns else "SNP"
        if snp_col_name not in clumped.columns:
            return pd.DataFrame(columns=df.columns)

        lead_rsids = set(clumped[snp_col_name].astype(str))

    return df[df[snp_col].isin(lead_rsids)].copy()


def find_proxies(
    missing_rsids: list[str],
    bfile: Path = LD_REF_PREFIX,
    ld_window_kb: int = 5000,
    r2_threshold: float = 0.8,
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
        ordered_targets = list(dict.fromkeys(str(rsid) for rsid in missing_rsids if str(rsid)))
        snplist.write_text("\n".join(ordered_targets))

        out_prefix = tmp / "ld"
        cmd = [
            _PLINK2,
            "--bfile",
            str(bfile),
            "--r2-unphased",
            "--ld-snp-list",
            str(snplist),
            "--ld-window-kb",
            str(ld_window_kb),
            "--ld-window-r2",
            str(r2_threshold),
            "--out",
            str(out_prefix),
        ]
        _run(cmd)

        vcor_file = out_prefix.with_suffix(".vcor")
        if not vcor_file.exists():
            return {}

        ld_df = pd.read_csv(vcor_file, sep="\t")
        if ld_df.empty or "ID_A" not in ld_df.columns or "ID_B" not in ld_df.columns:
            return {}

        r2_col = _pick_r2_col(list(ld_df.columns))
        if not r2_col:
            return {}

        target_set = set(ordered_targets)
        proxies: dict[str, tuple[str, float]] = {}

        for _, row in ld_df.iterrows():
            target = str(row["ID_A"])
            proxy = str(row["ID_B"])
            if target not in target_set or proxy == target:
                continue

            try:
                r2 = float(row[r2_col])
            except (TypeError, ValueError):
                continue
            if r2 < r2_threshold:
                continue

            current = proxies.get(target)
            if current is None or r2 > current[1] or (r2 == current[1] and proxy < current[0]):
                proxies[target] = (proxy, r2)

    return proxies


_MAJOR_MINOR_RE_TMPL = (
    r"{variant} alleles:\s+"
    r"MAJOR = (?:[A-Za-z0-9_]+ = )?([A-Za-z]+)\s+"
    r"MINOR = (?:[A-Za-z0-9_]+ = )?([A-Za-z]+)"
)
_PHASE_RE = re.compile(r"Major alleles are (in phase|out of phase) with each other\.", re.IGNORECASE)


def in_phase_allele_map(
    target_rsid: str,
    proxy_rsid: str,
    bfile: Path = LD_REF_PREFIX,
) -> dict[str, str] | None:
    """
    Return target->proxy in-phase allele mapping based on plink2 --ld output.

    Returns None if PLINK output is unavailable or cannot be parsed.
    """
    cmd = [
        _PLINK2,
        "--bfile",
        str(bfile),
        "--ld",
        target_rsid,
        proxy_rsid,
    ]
    try:
        result = _run(cmd)
    except Exception:
        return None

    text = f"{result.stdout}\n{result.stderr}"
    target_match = re.search(_MAJOR_MINOR_RE_TMPL.format(variant=re.escape(target_rsid)), text, re.MULTILINE)
    proxy_match = re.search(_MAJOR_MINOR_RE_TMPL.format(variant=re.escape(proxy_rsid)), text, re.MULTILINE)
    phase_matches = _PHASE_RE.findall(text)

    if not target_match or not proxy_match or not phase_matches:
        return None

    target_major, target_minor = target_match.group(1).upper(), target_match.group(2).upper()
    proxy_major, proxy_minor = proxy_match.group(1).upper(), proxy_match.group(2).upper()

    if phase_matches[-1].lower().startswith("in"):
        return {
            target_major: proxy_major,
            target_minor: proxy_minor,
        }

    return {
        target_major: proxy_minor,
        target_minor: proxy_major,
    }


def r_square_matrix(snp_list: list[str], bfile: Path = LD_REF_PREFIX) -> pd.DataFrame:
    """
    Compute pairwise r² matrix for snp_list via plink2 --r square.
    Returns a DataFrame (index and columns = rsIDs).
    """
    with tempfile.TemporaryDirectory(prefix="ld_matrix_") as tmpdir:
        tmp = Path(tmpdir)
        snpfile = tmp / "snps.txt"
        snpfile.write_text("\n".join(snp_list))

        out_prefix = tmp / "ld"
        cmd = [
            _PLINK2,
            "--bfile",
            str(bfile),
            "--extract",
            str(snpfile),
            "--r",
            "square",
            "--out",
            str(out_prefix),
        ]
        _run(cmd)

        matrix_file = out_prefix.with_suffix(".r.square")
        if not matrix_file.exists():
            raise FileNotFoundError(f"plink2 r-square output not found: {matrix_file}")

        mat = pd.read_csv(matrix_file, sep="\t", header=None)
        mat.index = snp_list
        mat.columns = snp_list
        return mat

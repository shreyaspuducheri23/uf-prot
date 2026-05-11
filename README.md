**Proteome-wide Mendelian randomization to identify causal plasma proteins for uterine fibroids**

Two-sample MR using cis-pQTLs from...

1. ARIC — SomaScan, 7K, European + African American
2. deCODE — SomaScan, 35K, Icelandic
3. Fenland — SomaScan, 10K, European
4. UKB-PPP — Olink, 34K, European

...as exposures against the largest available fibroid GWAS (Kim et al. 2025, GCST90461957, 74K cases / 466K controls) as the outcome. Followed by colocalization analysis to validate. 

Corollary analysis may include PheWAS-based safety screening and druggability analysis to prioritize druggable targets with a downstream goal of nominating candidates for drug-coated embolic agents in UAE.

First, we want to restrict analyses to the european-ancestry. This is for maximizing instrument strength and statistical power for the proof-of-concept. Planned extensions will include multi-ancestry and african ancestry. 

## Data

### Outcome


| Dataset                            | ID           | Build  | N    | Location                                                                      |
| ---------------------------------- | ------------ | ------ | ---- | ----------------------------------------------------------------------------- |
| Kim et al. 2025 fibroid GWAS (EUR) | GCST90461958 | GRCh38 | 434K | `data/raw/kim_fibroid_gwas/GCST90461958.h.tsv.gz` (harmonised, tabix-indexed) |


### Exposures (pQTL summary statistics)


| Cohort  | Platform     | N    | Format                                                              | Location                              |
| ------- | ------------ | ---- | ------------------------------------------------------------------- | ------------------------------------- |
| ARIC EA | SomaScan 7K  | ≈7K  | PLINK2 `.glm.linear`, one file per SeqId                            | `data/raw/ARIC/EA/` (4,657 proteins)  |
| ARIC AA | SomaScan 7K  | ≈7K  | Same                                                                | `data/raw/ARIC/AA/` (4,657 proteins)  |
| deCODE  | SomaScan 35K | ≈35K | Single compressed sentinel pQTL file                                | `data/raw/deCODE/`                    |
| UKB-PPP | Olink 34K    | ≈34K | Synapse `syn51365303`; 2,940 `.tar` files × ≈547 MB = ≈1.6 TB total | Stream-only (see below)               |
| Fenland | SomaScan 10K | ≈10K | Synapse `syn51824537`; ≈5K proteins, 2 files each                   | Stream-only (see below)               |


`data/raw/ARIC/seqid.txt` maps SeqId → UniProt → gene symbol → chromosome → TSS (hg19).

### LD reference

1000 Genomes EUR panel (502 samples, ≈9M variants, MAF > 0.01) at `data/ld_ref/` (symlink to MR_IA). PLINK bed/bim/fam format, consumed via PLINK2. Used for LD clumping and proxy SNP search. TODO: consider replacing with a UKB EUR subset (≈50k samples) for better LD estimates, particularly for uncommon variants — requires UKB data access to generate.

### Synapse streaming

UKB-PPP and Fenland cannot be downloaded in full (≈1.6 TB+ each; local disk insufficient). Use the streaming approach from MR_IA: download one file at a time, extract the relevant cis-region in memory, save matched rows, delete. Peak disk ≈550 MB per worker. Synapse credentials are read from `~/.synapseConfig`. Install client with `uv pip install synapseclient`.

## Methods

*Each script is standalone and tackles a thematic step, named as such. Scripts lie in `./scripts/` and processed data in `./processed_data/`.*

### Genetic instrument selection

Perform LD clumping using a clumping window of 1 Mb, r2 threshold of 0.001, and genome-wide significance level of 5 × 10-8 to identify independent variants from each of the four original studies. Variants within 500 kb of the transcription start site of the protein-coding gene were considered cis-pQTLs. We used ancestry-matched reference panels as described below (European: 1000G EUR, 502 samples; TODO: replace with UKB 50k EUR). Only variants with a MAF greater than 1% in these ancestry-matched panels were retained for further analysis.

### LD-clumping

LD clumping with a clumping window of 1 Mb, r2 threshold of 0.001, and genome-wide significance level of 5 × 10-8. LD is computed against ancestry-matched reference panels. For European ancestry, we currently use the 1000 Genomes EUR panel (502 samples; TODO: replace with UKB 50k EUR subset for improved LD resolution). TODO: define panels for other ancestries.

### Two-sample Mendelian randomization

We estimated the effect of circulating plasma protein levels on meta-analyzed leiomyoma (fibroid) outcomes using proteome-wide two-sample MR in a genetically stratified manner (TwoSampleMR). MR relies on three instrumental-variable assumptions: (1) Relevance: the genetic instrument is associated with the exposure; (2) Independence: the instrument is independent of any confounders of the exposure-outcome relationship; and (3) Exclusion restriction: the instrument influences the outcome only through its effect on the exposure (i.e., no alternative causal pathways or horizontal pleiotropy). We excluded proteins within the major histocompatibility complex (MHC) due to the region's complex linkage disequilibrium structure, which can confound genetic association signals and complicate interpretation of causal relationships.

For European analyses, we utilized proteomic GWAS summary data from ARIC, deCODE, Fenland, and UKB-PPP studies as exposures and assessed the effects of proteins in each cohort on European fibroid GWAS.  In this context, the term "protein" refers to the aptamer (SomaScan) or antibodies (Olink) targeting the respective protein.
TODO: define the cohorts for other ancestries, where ARIC will play a special role due to african ancestry representation.

Summary statistics from protein GWAS (exposures) and fibroid outcome GWAS were harmonized using the harmonise_data() function. When a selected instrument is unavailable in the outcome dataset, we conduct a proxy SNP search using ancestry-matched LD reference panels from the instrument selection procedure. Proxy variants are identified using PLINK2 with parameters --ld-window-kb=5000 and --ld-window-r2=0.8. We retain proxies with MAFs ≤ 0.42. Instrument strength is evaluated using F-statistics, with values above 10 considered indicative of robust instruments and lower likelihood of weak instrument bias.

Causal estimates were derived using the mr() function and represent the odds ratio per 1 standard deviation increase in genetically predicted protein level. For proteins with one instrument, we applied the Wald ratio; for those with two or more instruments, we used an inverse-variance weighted random-effects model. We controlled for multiple testing by using a Benjamini-Hochberg false discovery rate (FDR) threshold of 5%, consistent with prior studies. However, we performed the correction within each proteomics cohort and at the single-trait level. A protein reaching FDR < 0.05 in any single stratum is reported. Ensure MR analyses adhered to STROBE-MR reporting guidelines.

### Sensitivity analyses with alternative MR methods

To ensure the robustness of causal estimates, we applied multiple sensitivity checks, including heterogeneity testing, alternative MR methods (weighted median, weighted mode, MR-Egger), and Steiger directionality testing. For proteins with two or more instruments, we calculated a heterogeneity P value using Cochran's Q (Q_pval) and I^2. Associations with I^2 ≥ 0.5 and Q_pval < 0.05 were considered heterogeneous. For proteins with three or more instruments, we additionally applied weighted median, weighted mode, and MR-Egger methods. Consistency in effect direction across these methods was required. Directional pleiotropy was tested using mr_pleiotropy_test() (significant if P < 0.05). We also used Steiger filtering (directionality_test()) to exclude variants suggesting reverse causation.

### Colocalization

We used the colocalization method, SharePro, to test whether protein levels and GWAS outcomes share causal variants for associations passing sensitivity analyses above. Analyses were conducted within 1 Mb of the lead cis-pQTL using default priors. Colocalization was defined as a posterior probability (PP.H4) ≥ 0.8. As a sensitivity check, we additionally ran `coloc.abf` on the same regions.

### Druggability assessment

We evaluated the druggability of putatively causal proteins using the druggable genome by Finan et al., DrugBank, and Open Targets.

---

## Implementation

### Quick start

```bash
# Install Python + R dependencies, download liftover chain, clone SharePro
bash scripts/00_setup/install.sh

# Smoke-test on ARIC EA (local, fast) end-to-end before running full pipeline
uv run python scripts/02_cis_pqtl_extract/aric.py --limit 10
uv run python scripts/03_clump/clump.py --cohort ARIC_EA --limit 10
uv run python scripts/04_liftover/instruments_to_hg38.py --cohort ARIC_EA --limit 10
uv run python scripts/05_harmonise/harmonise.py --cohort ARIC_EA --limit 10
Rscript scripts/06_mr/run_mr.R --cohort ARIC_EA --limit 10
```

### Full pipeline

Each script is standalone, resumable (checkpointed), progress-instrumented, and logs to `logs/`.

```bash
# 1. Parse Kim GWAS metadata and validate tabix lookup
uv run python scripts/01_outcome_prep/prep_kim.py

# 2. Extract cis-pQTLs (±500 kb of TSS, p < 5×10⁻⁸, MAF > 1%, no MHC)
uv run python scripts/02_cis_pqtl_extract/aric.py
uv run python scripts/02_cis_pqtl_extract/decode.py          # HTTP, sequential
uv run python scripts/02_cis_pqtl_extract/ukbppp.py --workers 4   # Synapse streaming
uv run python scripts/02_cis_pqtl_extract/fenland.py              # Synapse streaming

# 3. LD clumping (1 Mb window, r² < 0.001, p < 5×10⁻⁸) vs 1000G EUR
uv run python scripts/03_clump/clump.py --cohort all

# 4. Lift instrument positions from hg19 → GRCh38 (deCODE: pass-through)
uv run python scripts/04_liftover/instruments_to_hg38.py --cohort all

# 5. Join instruments with Kim outcome; proxy SNP search for absent variants
uv run python scripts/05_harmonise/harmonise.py --cohort all

# 6. Two-sample MR (Wald ratio / IVW-MRE) + BH-FDR within each cohort
Rscript scripts/06_mr/run_mr.R

# 7. Sensitivity analyses (Q, I², weighted median/mode, MR-Egger, Steiger)
Rscript scripts/07_sensitivity/run_sensitivity.R

# 8. Colocalization on FDR+sensitivity-passing proteins
uv run python scripts/08_coloc/extract_regions.py --cohort all   # ±1 Mb regions
uv run python scripts/08_coloc/sharepro.py --cohort all          # primary
Rscript scripts/08_coloc/coloc_abf.R                             # sensitivity

# 9. Assemble tiered final results table
uv run python scripts/09_assemble/assemble.py
```

All Python scripts accept `--limit N` to cap the number of proteins (for testing). Scripts that have already completed proteins skip them on re-run — Ctrl-C safe.

### Genome builds

Each cohort's extraction and clumping operate in the cohort's native build. Step 4 lifts surviving instrument SNPs to GRCh38 to align with the Kim outcome (GRCh38). The 1000G EUR LD reference (GRCh37) is used only for clumping and proxy search, both of which operate in native coords before liftover.

| Source          | Build   | Lifted in step |
|-----------------|---------|----------------|
| ARIC EA         | GRCh37  | 04             |
| deCODE          | GRCh38  | pass-through   |
| UKB-PPP         | GRCh37  | 04             |
| Fenland         | GRCh37  | 04             |
| Kim (outcome)   | GRCh38  | — (reference)  |
| 1000G EUR LD ref | GRCh37 | clump/proxy only |

### Code layout

```
scripts/
  lib/                      # shared Python utilities (imported by all step scripts)
    paths.py                # canonical paths — single source of truth
    filters.py              # MAF, MHC exclusion, palindrome, GW significance
    cis.py                  # ±500 kb / ±1 Mb window bounds; TSS lookup (Ensembl REST)
    cis_extract.py          # cohort-agnostic extraction loop (shared by all 4 cohorts)
    plink.py                # PLINK2 subprocess wrappers: clump, proxies, LD matrix
    liftover.py             # hg19 ↔ hg38 coordinate conversion (pyliftover)
    outcome.py              # Kim GWAS tabix lookup (OutcomeLookup class)
    synapse_stream.py       # tar/gz streaming helpers for UKB-PPP and Fenland
    decode_stream.py        # HTTP download helpers for deCODE signed URLs
    sumstats_io.py          # read/write normalised TSVs with canonical dtypes
    fstat.py                # F-statistic computation
    fdr.py                  # Benjamini-Hochberg FDR
    checkpoint.py           # per-unit resumable state ledger (JSON)
    logging.py              # ISO-timestamped logger + run manifest
    progress.py             # tqdm wrapper with consistent formatting
    schema.py               # canonical column names; ProteinMeta dataclass

  rlib/                     # shared R utilities (sourced by R step scripts)
    harmonise.R             # TwoSampleMR::harmonise_data wrapper
    mr_methods.R            # Wald / IVW-MRE / sensitivity analysis drivers
    coloc_abf.R             # coloc.abf wrapper
    logging.R               # futile.logger with ISO timestamps + file sink
    checkpoint.R            # RDS-based per-unit checkpointing
    progress.R              # progressr / pbmcapply wrappers

  00_setup/                 # install.sh, install_packages.R, install_sharepro.sh
  01_outcome_prep/          # prep_kim.py
  02_cis_pqtl_extract/      # aric.py, decode.py, ukbppp.py, fenland.py
  03_clump/                 # clump.py
  04_liftover/              # instruments_to_hg38.py
  05_harmonise/             # harmonise.py
  06_mr/                    # run_mr.R
  07_sensitivity/           # run_sensitivity.R
  08_coloc/                 # extract_regions.py, sharepro.py, coloc_abf.R
  09_assemble/              # assemble.py

processed_data/
  {ARIC_EA,deCODE,UKB_PPP,Fenland}/
    cis_sumstats/           # per-protein filtered cis-pQTL TSVs (step 02)
    instruments/            # post-clumping independent instruments (step 03)
    instruments_hg38/       # after GRCh38 liftover (step 04)
    harmonised/             # after harmonisation with Kim (step 05)
    mr_results.tsv          # MR estimates + BH-FDR (step 06)
    sensitivity.tsv         # sensitivity results + passes_sensitivity flag (step 07)
  outcome/                  # Kim metadata JSON
  coloc/                    # SharePro + coloc.abf results by cohort
  mr_all_cohorts.tsv        # combined MR results across cohorts
  final_results.tsv         # tiered output (step 09)

tools/
  SharePro_coloc/           # cloned from github.com/zhwm/SharePro_coloc at setup

data/
  raw/                      # input summary statistics (see Data section above)
  ld_ref/                   # 1000G EUR PLINK binary (symlink)
  ref/                      # hg19ToHg38.over.chain.gz (downloaded at setup)

logs/
  _manifest.tsv             # one-line record per script run (script, args, timing, n_units)
  <step>_<timestamp>.log    # full timestamped log per run
```

### Results tiering

`processed_data/final_results.tsv` assigns each protein-cohort association to a tier:

| Tier               | Criteria                                                                       |
|--------------------|--------------------------------------------------------------------------------|
| Tier 1 (replicated)| FDR < 0.05 + passes sensitivity + SharePro PP.H4 ≥ 0.8 + coloc.abf agrees   |
| Tier 1             | FDR < 0.05 + passes sensitivity + SharePro PP.H4 ≥ 0.8                       |
| Tier 2             | FDR < 0.05 + passes sensitivity (colocalization inconclusive)                 |
| Tier 2 (no sens)   | FDR < 0.05, sensitivity not computable (single-SNP Wald ratio instruments)    |

### Tests

```bash
uv run pytest tests/ -v
```

85 tests cover all shared library modules (filters, FDR, F-stat, checkpointing, sumstats I/O, cis-window logic, ARIC file parsing, outcome tabix lookup, Synapse streaming, deCODE streaming, extraction pipeline). Tests requiring absent files (liftover chain, Synapse credentials) are automatically skipped.

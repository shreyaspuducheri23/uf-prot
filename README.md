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
| Kim et al. 2025 fibroid GWAS (EUR) | GCST90461958 | GRCh37 | 434K | `data/raw/kim_fibroid_gwas/GCST90461958.h.tsv.gz` (harmonised, tabix-indexed) |


### Exposures (pQTL summary statistics)


| Cohort  | Platform     | N    | Format                                                              | Location                              |
| ------- | ------------ | ---- | ------------------------------------------------------------------- | ------------------------------------- |
| ARIC EA | SomaScan 7K  | ≈7K  | PLINK2 `.glm.linear`, one file per SeqId                            | `data/raw/ARIC/EA/` (4,657 proteins)  |
| ARIC AA | SomaScan 7K  | ≈7K  | Same                                                                | `data/raw/ARIC/AA/` (4,657 proteins)  |
| deCODE  | SomaScan 35K | ≈35K | Single compressed sentinel pQTL file                                | `data/raw/deCODE/` (symlink to MR_IA) |
| UKB-PPP | Olink 34K    | ≈34K | Synapse `syn51365303`; 2,940 `.tar` files × ≈547 MB = ≈1.6 TB total | Stream-only (see below)               |
| Fenland | SomaScan 10K | ≈10K | Synapse `syn51824537`; ≈5K proteins, 2 files each                   | Stream-only (see below)               |


`data/raw/ARIC/seqid.txt` maps SeqId → UniProt → gene symbol → chromosome → TSS (hg19).

### LD reference

1000 Genomes EUR panel (502 samples, ≈9M variants, MAF > 0.01) at `data/ld_ref/` (symlink to MR_IA). PLINK1 binary format. Used for LD clumping and HEIDI. TODO: consider replacing with a UKB EUR subset (≈50k samples) for better LD estimates, particularly for uncommon variants — requires UKB data access to generate.

### Synapse streaming

UKB-PPP and Fenland cannot be downloaded in full (≈1.6 TB+ each; local disk insufficient). Use the streaming approach from MR_IA: download one file at a time, extract the relevant cis-region in memory, save matched rows, delete. Peak disk ≈550 MB per worker. Synapse credentials are read from `~/.synapseConfig`. Install client with `uv pip install synapseclient`.

## Methods

*each script should standalone and should tackle a thematic step, and should be named as such. scripts should lie in ./code/ and processed data in ./processed_data/*

### Genetic instrument selection

perform LD clumping using a clumping window of 1 Mb, r2 threshold of 0.001, and genome-wide significance level of 5 × 10-8 to identify independent variants from each of the four original studies. Variants within 500 kb of the transcription start site of the protein-coding gene were considered cis-pQTLs. We used ancestry-matched reference panels as described below (European: 1000G EUR, 502 samples; TODO: replace with UKB 50k EUR). Only variants with a MAF greater than 1% in these ancestry-matched panels were retained for further analysis.

### LD-clumping

LD clumping with a clumping window of 1 Mb, r2 threshold of 0.001, and genome-wide significance level of 5 × 10-8. LD is computed against ancestry-matched reference panels. For European ancestry, we currently use the 1000 Genomes EUR panel (502 samples; TODO: replace with UKB 50k EUR subset for improved LD resolution). TODO: define panels for other ancestries.

### Two-sample Mendelian randomization

We estimated the effect of circulating plasma protein levels on meta-analyzed leiomyoma (fibroid) outcomes using proteome-wide two-sample MR in a genetically stratified manner (TwoSampleMR). MR relies on three instrumental-variable assumptions: (1) Relevance: the genetic instrument is associated with the exposure; (2) Independence: the instrument is independent of any confounders of the exposure-outcome relationship; and (3) Exclusion restriction: the instrument influences the outcome only through its effect on the exposure (i.e., no alternative causal pathways or horizontal pleiotropy). We excluded proteins within the major histocompatibility complex (MHC) due to the region's complex linkage disequilibrium structure, which can confound genetic association signals and complicate interpretation of causal relationships.

For European analyses, we utilized proteomic GWAS summary data from ARIC, deCODE, Fenland, and UKB-PPP studies as exposures and assessed the effects of proteins in each cohort on European fibroid GWAS.  In this context, the term "protein" refers to the aptamer (SomaScan) or antibodies (Olink) targeting the respective protein.
TODO: define the cohorts for other ancestries, where ARIC will play a special role due to african ancestry representation.

Summary statistics from protein GWAS (exposures) and fibroid outcome GWAS were harmonized using the harmonise_data() function. When a selected instrument is unavailable in the outcome dataset, we conduct a proxy SNP search using ancestry-matched LD reference panels from the instrument selection procedure. Proxy variants are identified using PLINK v.1.944 with parameters --ld-window=5000, --ld-window-kb=5000, --ld-window-r2=0.8. We retain proxies with MAFs ≤ 0.42. Instrument strength is evaluated using F-statistics, with values above 10 considered indicative of robust instruments and lower likelihood of weak instrument bias.

Causal estimates were derived using the mr() function and represent the odds ratio per 1 standard deviation increase in genetically predicted protein level. For proteins with one instrument, we applied the Wald ratio; for those with two or more instruments, we used an inverse-variance weighted random-effects model. We controlled for multiple testing by using a Benjamini-Hochberg false discovery rate (FDR) threshold of 5%, consistent with prior studies. However, we performed the correction within each proteomics cohort and at the single-trait level. A protein reaching FDR < 0.05 in any single stratum is reported. Ensure MR analyses adhered to STROBE-MR reporting guidelines.

### Sensitivity analyses with alternative MR methods

To ensure the robustness of causal estimates, we applied multiple sensitivity checks, including heterogeneity testing, alternative MR methods (weighted median, weighted mode, MR-Egger), and Steiger directionality testing. For proteins with two or more instruments, we calculated a heterogeneity P value using Cochran’s Q (Q_pval) and I^2. Associations with I^2 ≥ 0.5 and Q_pval < 0.05 were considered heterogeneous. For proteins with three or more instruments, we additionally applied weighted median, weighted mode, and MR-Egger methods. Consistency in effect direction across these methods was required. Directional pleiotropy was tested using mr_pleiotropy_test() (significant if P < 0.05). We also used Steiger filtering (directionality_test()) to exclude variants suggesting reverse causation.

### Colocalization

We used the colocalization method, SharePro, to test whether protein levels and GWAS outcomes share causal variants for associations passing sensitivity analyses above. Analyses were conducted within 1 Mb of the lead cis-pQTL using default priors. Colocalization was defined as a posterior probability (PP.H4) ≥ 0.8.

### Druggability assessment

We evaluated the druggability of putatively causal proteins using the druggable genome by Finan et al., DrugBank, and Open Targets.
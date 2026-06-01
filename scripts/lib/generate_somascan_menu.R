#!/usr/bin/env Rscript
# Generate data/raw/somascan/analyte_menu.tsv from the SomaScan.db Bioconductor package.
# Run once from the repo root:
#   Rscript scripts/lib/generate_somascan_menu.R
# Requires: BiocManager::install("SomaScan.db")

library(SomaScan.db)

all_keys <- keys(SomaScan.db)  # "14157-21" hyphenated format

df <- select(SomaScan.db,
             keys    = all_keys,
             columns = c("PROBEID", "SYMBOL"),
             keytype = "PROBEID")

# Convert from hyphen (14157-21) to underscore (14157_21) to match Fenland/deCODE filenames
df$seqid_key   <- gsub("-", "_", df$PROBEID)
df$gene_symbol <- df$SYMBOL

out <- df[!is.na(df$SYMBOL), c("seqid_key", "gene_symbol")]

dest <- "data/raw/somascan/analyte_menu.tsv"
dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
write.table(out, dest, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote %d analytes to %s\n", nrow(out), dest))

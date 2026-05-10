# Install all required R packages
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes",    quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("BiocManager",quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c(
  "data.table", "dplyr",
  "coloc",
  "progress",
  "pbmcapply",
  "futile.logger"
)
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  remotes::install_github("MRCIEU/TwoSampleMR")
}
if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  remotes::install_github("rondolab/MR-PRESSO")
}

cat("All R packages installed.\n")

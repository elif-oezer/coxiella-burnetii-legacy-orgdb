# 00_setup_environment.R
#
# Install the R/Bioconductor packages required by the OrgDb build workflow.
#
# Run this script manually when setting up the project environment.

#Install BioManager if needed:

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

#Bioconductor dependencies
BiocManager::install(c(
  "AnnotationForge",
  "AnnotationDbi",
  "rtracklayer",
  "Rsamtools"
))

#CRAN dependencies
install.packages(c(
  "dplyr",
  "here"
))



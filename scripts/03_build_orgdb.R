# ==============================================================================
# 03_build_orgdb.R
#
# Build a custom OrgDb package for Coxiella burnetii using AnnotationForge.
#
# Inputs:
#   - output/gene_info.rds
#   - output/go_final.rds
#
# Output:
#   - generated OrgDb package source directory inside output/
#
# Primary identifier:
#   GID = historical locus tag (old_locus_tag)
#
# Additional identifier:
#   LOCUS_TAG = current locus tag
#
# ==============================================================================


# ---- Packages ---------------------------------------------------------------

library(AnnotationForge)
library(AnnotationDbi)
library(dplyr)

source(here::here("scripts", "00_config.R"))


# ---- Check required prepared inputs -----------------------------------------

if (!file.exists(path_gene_info_rds)) {
  stop(
    "Prepared gene annotation file not found: ",
    path_gene_info_rds,
    "\nRun scripts/01_prepare_annotation_data.R first."
  )
}

if (!file.exists(path_go_final_rds)) {
  stop(
    "Prepared GO annotation file not found: ",
    path_go_final_rds,
    "\nRun scripts/01_prepare_annotation_data.R first."
  )
}


# ---- Load prepared annotation tables ----------------------------------------

message("Loading prepared annotation tables...")

gene_info <- readRDS(path_gene_info_rds)

go_final <- readRDS(path_go_final_rds)


message(
  "Loaded ",
  nrow(gene_info),
  " genes and ",
  nrow(go_final),
  " GO associations."
)


# ---- Validate required columns ----------------------------------------------

required_gene_columns <- c(
  "GID",
  "SYMBOL",
  "GENENAME",
  "LOCUS_TAG"
)

missing_gene_columns <- setdiff(
  required_gene_columns,
  colnames(gene_info)
)

if (length(missing_gene_columns) > 0) {
  stop(
    "Required gene_info columns are missing: ",
    paste(missing_gene_columns, collapse = ", ")
  )
}


required_go_columns <- c(
  "GID",
  "GO",
  "EVIDENCE"
)

missing_go_columns <- setdiff(
  required_go_columns,
  colnames(go_final)
)

if (length(missing_go_columns) > 0) {
  stop(
    "Required GO columns are missing: ",
    paste(missing_go_columns, collapse = ", ")
  )
}


# ---- Validate gene identifiers ----------------------------------------------

if (anyNA(gene_info$GID)) {
  stop("gene_info contains missing GID values.")
}

if (any(gene_info$GID == "")) {
  stop("gene_info contains empty GID values.")
}

if (anyDuplicated(gene_info$GID)) {
  stop("gene_info contains duplicated GID values.")
}


# ---- Validate GO identifiers -------------------------------------------------

if (!all(go_final$GID %in% gene_info$GID)) {

  missing_ids <- setdiff(
    unique(go_final$GID),
    unique(gene_info$GID)
  )

  stop(
    "GO table contains GIDs not present in gene_info: ",
    paste(head(missing_ids, 20), collapse = ", ")
  )
}


# ---- Prepare tables for AnnotationForge --------------------------------------
#
# makeOrgPackage() expects:
#
#   gene_info:
#       GID + arbitrary annotation columns
#
#   go:
#       GID
#       GO
#       EVIDENCE
#
# Additional columns such as LOCUS_TAG become accessible through the OrgDb.

gene_info_build <- gene_info %>%

  select(
    GID,
    SYMBOL,
    GENENAME,
    LOCUS_TAG
  ) %>%

  distinct()


go_build <- go_final %>%

  select(
    GID,
    GO,
    EVIDENCE
  ) %>%

  distinct()


# ---- Package name ------------------------------------------------------------

pkg_name <- paste0(
  "org.",
  substr(orgdb_genus, 1, 1),
  orgdb_species,
  ".eg.db"
)

pkg_path <- file.path(
  dir_orgdb_build,
  pkg_name
)


# ---- Remove previous generated package --------------------------------------

if (dir.exists(pkg_path)) {

  message(
    "Removing existing generated package directory: ",
    pkg_path
  )

  unlink(
    pkg_path,
    recursive = TRUE,
    force = TRUE
  )
}


# ---- Build OrgDb package -----------------------------------------------------

message("Building OrgDb package: ", pkg_name)
if (!grepl(
  "^[^<>]+ <[^<>[:space:]]+@[^<>[:space:]]+\\.[^<>[:space:]]+>$",
  orgdb_maintainer
)) {
  stop(
    "Invalid orgdb_maintainer format: ",
    orgdb_maintainer,
    "\nExpected: Name <email@example.com>"
  )
}

message("Maintainer passed to AnnotationForge: ", orgdb_maintainer)
AnnotationForge::makeOrgPackage(

  gene_info = gene_info_build,

  go = go_build,

  version = orgdb_version,

  maintainer = orgdb_maintainer,

  author = orgdb_author,

  outputDir = dir_orgdb_build,

  tax_id = orgdb_tax_id,

  genus = orgdb_genus,

  species = orgdb_species,

  goTable = "go"
)


# ---- Confirm package creation ------------------------------------------------

if (!dir.exists(pkg_path)) {
  stop(
    "OrgDb build finished but package directory was not created: ",
    pkg_path
  )
}


# ---- Summary ----------------------------------------------------------------

message(
  "\nOrgDb package build complete.\n",
  "\n",
  "Package name: ",
  pkg_name,
  "\n",
  "Package directory: ",
  pkg_path,
  "\n",
  "Genes included: ",
  nrow(gene_info_build),
  "\n",
  "GO associations included: ",
  nrow(go_build),
  "\n",
  "Package version: ",
  orgdb_version
)

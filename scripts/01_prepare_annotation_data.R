# ==============================================================================
# 01_prepare_annotation_data.R
#
# Prepare annotation tables required by AnnotationForge::makeOrgPackage().
#
# Inputs:
#   - NCBI genome annotation GFF
#   - NCBI/RefSeq Gene Ontology annotation table
#
# Outputs:
#   - output/gene_info.rds
#   - output/go_final.rds
#   - output/coxiella_go_clean.tsv
#
# Identifier strategy:
#
#   GFF old_locus_tag  -> GID
#   GO  Locus_tag      -> GID
#
# The current GFF locus_tag is retained as an additional annotation field.
#
# Validation is performed separately in 02_qc_checks.R.
# ==============================================================================


# ---- Packages ---------------------------------------------------------------

library(rtracklayer)
library(dplyr)

source(here::here("scripts", "00_config.R"))


# ---- Check required input files ---------------------------------------------

if (!file.exists(path_gff)) {
  stop("Genome annotation file not found: ", path_gff)
}

if (!file.exists(path_go_raw)) {
  stop("GO annotation file not found: ", path_go_raw)
}


# ---- Helper function ---------------------------------------------------------
#
# Return the first non-missing, non-empty value in a vector.
# This is safer than first(na.omit(x)) when all values are missing.

first_nonempty <- function(x) {

  x <- as.character(x)

  x <- x[
    !is.na(x) &
      x != ""
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[[1]]
}


# ---- Import genome annotation -----------------------------------------------

message("Importing genome annotation...")

gff <- rtracklayer::import(path_gff)

gff_df <- as.data.frame(gff)


# ---- Check required GFF columns ---------------------------------------------

required_gff_columns <- c(
  "type",
  "old_locus_tag",
  "locus_tag"
)

missing_gff_columns <- setdiff(
  required_gff_columns,
  colnames(gff_df)
)

if (length(missing_gff_columns) > 0) {
  stop(
    "Required GFF columns are missing: ",
    paste(missing_gff_columns, collapse = ", ")
  )
}


# Some metadata fields may not occur in every GFF.
# Create them as NA if absent so downstream code remains robust.

optional_gff_columns <- c(
  "Name",
  "gene",
  "product",
  "gene_biotype"
)

for (col in optional_gff_columns) {

  if (!col %in% colnames(gff_df)) {
    gff_df[[col]] <- NA_character_
  }
}


# ---- Prepare gene annotation table ------------------------------------------
#
# The GO annotation uses the historical locus-tag namespace.
# Therefore old_locus_tag is used as the primary package gene identifier (GID).
#
# Current locus_tag values are retained as metadata.
#
# Both gene and CDS features are used because useful annotation fields may
# occur on either record type.

message("Preparing gene information table...")


# ---- Gene features -----------------------------------------------------------

genes_tbl <- gff_df %>%

  filter(type == "gene") %>%

  transmute(

    GID = as.character(old_locus_tag),

    SYMBOL = coalesce(
      as.character(gene),
      as.character(Name),
      as.character(old_locus_tag)
    ),

    GENENAME = coalesce(
      as.character(product),
      as.character(gene_biotype),
      as.character(Name),
      as.character(old_locus_tag)
    ),

    LOCUS_TAG = as.character(locus_tag),

    OLD_LOCUS_TAG = as.character(old_locus_tag)
  ) %>%

  filter(
    !is.na(GID),
    GID != ""
  ) %>%

  distinct()


# ---- CDS features ------------------------------------------------------------

cds_tbl <- gff_df %>%

  filter(type == "CDS") %>%

  transmute(

    GID = as.character(old_locus_tag),

    SYMBOL = coalesce(
      as.character(gene),
      as.character(Name),
      as.character(old_locus_tag)
    ),

    GENENAME = coalesce(
      as.character(product),
      as.character(Name),
      as.character(old_locus_tag)
    ),

    LOCUS_TAG = as.character(locus_tag),

    OLD_LOCUS_TAG = as.character(old_locus_tag)
  ) %>%

  filter(
    !is.na(GID),
    GID != ""
  ) %>%

  distinct()


# ---- Collapse to one row per GID --------------------------------------------

gene_info <- bind_rows(
  genes_tbl,
  cds_tbl
) %>%

  group_by(GID) %>%

  summarise(

    SYMBOL = first_nonempty(SYMBOL),

    GENENAME = first_nonempty(GENENAME),

    LOCUS_TAG = first_nonempty(LOCUS_TAG),

    OLD_LOCUS_TAG = first_nonempty(OLD_LOCUS_TAG),

    .groups = "drop"
  )


message(
  "Gene information table contains ",
  nrow(gene_info),
  " unique GIDs."
)


# ---- Import original RefSeq GO annotation -----------------------------------

message("Importing RefSeq GO annotations...")

go_raw <- read.delim(
  path_go_raw,
  header = TRUE,
  sep = "\t",
  quote = "\"",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

message(
  "GO annotation imported: ",
  nrow(go_raw),
  " rows and ",
  ncol(go_raw),
  " columns."
)


# ---- Check required GO columns ----------------------------------------------
#
# IMPORTANT:
#
# The GO file still uses the column name "Locus_tag".
#
# We are NOT expecting "old_locus_tag" in the GO file.
# Its Locus_tag values correspond to the GFF old_locus_tag namespace.

required_go_columns <- c(
  "Locus_tag",
  "GO_ID",
  "Evidence_Code"
)

missing_go_columns <- setdiff(
  required_go_columns,
  colnames(go_raw)
)

if (length(missing_go_columns) > 0) {
  stop(
    "Required GO columns are missing: ",
    paste(missing_go_columns, collapse = ", ")
  )
}


# ---- Prepare AnnotationForge GO table ---------------------------------------
#
# Original GO columns:
#
#   Locus_tag      -> GID
#   GO_ID          -> GO
#   Evidence_Code  -> EVIDENCE

go <- go_raw %>%

  transmute(

    GID = as.character(Locus_tag),

    GO = as.character(GO_ID),

    EVIDENCE = as.character(Evidence_Code)
  ) %>%

  filter(

    !is.na(GID),
    GID != "",

    !is.na(GO),
    GO != "",
    grepl("^GO:[0-9]{7}$", GO),

    !is.na(EVIDENCE),
    EVIDENCE != ""
  ) %>%

  distinct()


# ---- Save clean GO mapping ---------------------------------------------------

write.table(
  go,
  path_go_clean_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ---- Match GO annotation to genome annotation -------------------------------

go_final <- go %>%

  filter(
    GID %in% gene_info$GID
  ) %>%

  distinct()


# ---- Save prepared annotation objects ---------------------------------------

saveRDS(
  gene_info,
  path_gene_info_rds
)

saveRDS(
  go_final,
  path_go_final_rds
)


# ---- Summary ----------------------------------------------------------------

n_gene_ids <- n_distinct(gene_info$GID)

n_go_ids <- n_distinct(go$GID)

n_matching_go_ids <- length(
  intersect(
    unique(go$GID),
    unique(gene_info$GID)
  )
)

n_unmatched_go_ids <- length(
  setdiff(
    unique(go$GID),
    unique(gene_info$GID)
  )
)


message(
  "\nAnnotation preparation complete.\n",
  "\n",
  "Genes in gene_info: ",
  n_gene_ids,
  "\n",
  "GO associations in source annotation: ",
  nrow(go),
  "\n",
  "Unique GO-annotated GIDs: ",
  n_go_ids,
  "\n",
  "GO GIDs matching genome annotation: ",
  n_matching_go_ids,
  "\n",
  "GO GIDs not found in genome annotation: ",
  n_unmatched_go_ids,
  "\n",
  "GO associations retained after GID matching: ",
  nrow(go_final),
  "\n",
  "Unique GO-annotated genes retained: ",
  n_distinct(go_final$GID)
)

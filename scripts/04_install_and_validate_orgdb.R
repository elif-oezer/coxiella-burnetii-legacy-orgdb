# ==============================================================================
# 04_install_and_validate_orgdb.R
#
# Install and validate the generated Coxiella burnetii OrgDb package.
#
# Inputs:
#   - generated OrgDb package directory from 03_build_orgdb.R
#
# Validates:
#   - package installation
#   - package loading
#   - available keytypes and columns
#   - GID mappings
#   - current LOCUS_TAG mappings
#   - GO mappings
#
# ==============================================================================


# ---- Packages ---------------------------------------------------------------

library(AnnotationDbi)

source(here::here("scripts", "00_config.R"))


# ---- Package name and path ---------------------------------------------------

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


# ---- Check package source directory ------------------------------------------

if (!dir.exists(pkg_path)) {
  stop(
    "Generated OrgDb package directory not found: ",
    pkg_path,
    "\nRun scripts/03_build_orgdb.R first."
  )
}


# ---- Install generated package ----------------------------------------------

message("Installing OrgDb package from: ", pkg_path)

install.packages(
  pkg_path,
  repos = NULL,
  type = "source"
)


# ---- Load package ------------------------------------------------------------

message("Loading installed package: ", pkg_name)

suppressPackageStartupMessages(
  library(
    pkg_name,
    character.only = TRUE
  )
)


# ---- Retrieve OrgDb object ---------------------------------------------------

orgdb <- get(pkg_name)


# ---- Basic metadata ----------------------------------------------------------

message("\nPACKAGE METADATA")

print(orgdb)


# ---- Available keytypes ------------------------------------------------------

available_keytypes <- AnnotationDbi::keytypes(orgdb)

message("\nAVAILABLE KEYTYPES")

print(available_keytypes)


# ---- Available columns -------------------------------------------------------

available_columns <- AnnotationDbi::columns(orgdb)

message("\nAVAILABLE COLUMNS")

print(available_columns)


# ---- Number of primary keys --------------------------------------------------

all_gids <- AnnotationDbi::keys(
  orgdb,
  keytype = "GID"
)

message(
  "\nNumber of GID keys in OrgDb: ",
  length(all_gids)
)


# ---- Test GID lookup ---------------------------------------------------------

test_gids <- head(
  all_gids,
  5
)

message("\nTEST: GID -> annotation")

gid_test <- AnnotationDbi::select(
  orgdb,
  keys = test_gids,
  keytype = "GID",
  columns = intersect(
    c(
      "GID",
      "SYMBOL",
      "GENENAME",
      "LOCUS_TAG"
    ),
    available_columns
  )
)

print(gid_test)


# ---- Test LOCUS_TAG lookup ---------------------------------------------------

if ("LOCUS_TAG" %in% available_keytypes) {

  locus_tags <- AnnotationDbi::keys(
    orgdb,
    keytype = "LOCUS_TAG"
  )

  test_locus_tags <- head(
    locus_tags,
    5
  )

  message("\nTEST: LOCUS_TAG -> GID")

  locus_test <- AnnotationDbi::select(
    orgdb,
    keys = test_locus_tags,
    keytype = "LOCUS_TAG",
    columns = intersect(
      c(
        "GID",
        "SYMBOL",
        "GENENAME",
        "LOCUS_TAG"
      ),
      available_columns
    )
  )

  print(locus_test)

} else {

  warning(
    "LOCUS_TAG is not available as a keytype."
  )
}


# ---- Test GO mappings --------------------------------------------------------

if ("GO" %in% available_columns) {

  message("\nTEST: GID -> GO")

  go_test <- AnnotationDbi::select(
    orgdb,
    keys = test_gids,
    keytype = "GID",
    columns = intersect(
      c(
        "GO",
        "EVIDENCE",
        "ONTOLOGY"
      ),
      available_columns
    )
  )

  print(go_test)

} else {

  warning(
    "GO is not available as an OrgDb column."
  )
}


# ---- Count genes with GO annotations ----------------------------------------

if ("GO" %in% available_columns) {

  go_map <- AnnotationDbi::select(
    orgdb,
    keys = all_gids,
    keytype = "GID",
    columns = "GO"
  )

  go_map <- go_map[
    !is.na(go_map$GO),
    ,
    drop = FALSE
  ]

  message(
    "\nGenes with >=1 GO annotation in installed OrgDb: ",
    length(unique(go_map$GID))
  )

  message(
    "Total GID-GO associations in installed OrgDb: ",
    nrow(
      unique(
        go_map[, c("GID", "GO")]
      )
    )
  )
}


# ---- Simple consistency checks ----------------------------------------------

if (length(all_gids) == 0) {
  stop(
    "Validation failed: no GID keys were found."
  )
}


if (!"GID" %in% available_keytypes) {
  stop(
    "Validation failed: GID is not available as a keytype."
  )
}


if (!"LOCUS_TAG" %in% available_columns) {
  warning(
    "LOCUS_TAG is not available as an annotation column."
  )
}


if (!"GO" %in% available_columns) {
  stop(
    "Validation failed: GO annotations are unavailable."
  )
}


# ---- Validation summary ------------------------------------------------------

message(
  "\n============================================================\n",
  "ORGDB VALIDATION COMPLETE\n",
  "============================================================\n",
  "Package: ",
  pkg_name,
  "\n",
  "Version: ",
  as.character(
    utils::packageVersion(pkg_name)
  ),
  "\n",
  "GID keys: ",
  length(all_gids),
  "\n",
  "LOCUS_TAG keytype available: ",
  "LOCUS_TAG" %in% available_keytypes,
  "\n",
  "GO column available: ",
  "GO" %in% available_columns,
  "\n",
  "============================================================"
)

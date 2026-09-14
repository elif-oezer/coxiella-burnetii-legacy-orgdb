# ==============================================================================
# 05_clusterprofiler_smoke_test.R
#
# Functional smoke test of the custom Coxiella burnetii OrgDb package.
#
# Purpose:
#   Verify that the generated OrgDb can be used by clusterProfiler for
#   GO over-representation analysis (ORA).
#
# This is a technical validation, NOT a biological analysis.
# ==============================================================================


# ---- Packages ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clusterProfiler)
})

source(here::here("scripts", "00_config.R"))


# ---- OrgDb package -----------------------------------------------------------

pkg_name <- paste0(
  "org.",
  substr(orgdb_genus, 1, 1),
  orgdb_species,
  ".eg.db"
)


if (!requireNamespace(pkg_name, quietly = TRUE)) {
  stop(
    "OrgDb package is not installed: ",
    pkg_name,
    "\nRun scripts/04_install_and_validate_orgdb.R first."
  )
}


suppressPackageStartupMessages(
  library(
    pkg_name,
    character.only = TRUE
  )
)


OrgDb <- get(pkg_name)


message("Using OrgDb package: ", pkg_name)
message("Version: ", as.character(packageVersion(pkg_name)))


# ---- Confirm required keytypes/columns --------------------------------------

required_keytypes <- c("GID")

missing_keytypes <- setdiff(
  required_keytypes,
  AnnotationDbi::keytypes(OrgDb)
)

if (length(missing_keytypes) > 0) {
  stop(
    "Required OrgDb keytypes missing: ",
    paste(missing_keytypes, collapse = ", ")
  )
}


required_columns <- c("GO", "ONTOLOGY")

missing_columns <- setdiff(
  required_columns,
  AnnotationDbi::columns(OrgDb)
)

if (length(missing_columns) > 0) {
  stop(
    "Required OrgDb columns missing: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ---- Retrieve GO mappings ---------------------------------------------------

message("\nRetrieving GO annotations from OrgDb...")

all_gids <- AnnotationDbi::keys(
  OrgDb,
  keytype = "GID"
)


go_map <- AnnotationDbi::select(
  OrgDb,
  keys = all_gids,
  keytype = "GID",
  columns = c(
    "GO",
    "ONTOLOGY"
  )
)


go_map <- go_map[
  !is.na(go_map$GO) &
    !is.na(go_map$ONTOLOGY),
  ,
  drop = FALSE
]


message("Total GIDs in OrgDb: ", length(all_gids))

message(
  "GIDs with GO annotation: ",
  length(unique(go_map$GID))
)

message(
  "Unique GO terms represented: ",
  length(unique(go_map$GO))
)


# ---- Construct reproducible technical test gene set -------------------------
#
# We deliberately create a GO-biased gene set so that enrichment is likely
# to be detected.
#
# This is ONLY a functional test of the OrgDb/clusterProfiler interface.

go_counts <- table(go_map$GO)

candidate_terms <- names(
  sort(
    go_counts,
    decreasing = TRUE
  )
)


test_go <- NULL
test_genes <- NULL


for (term in candidate_terms) {

  genes_for_term <- unique(
    go_map$GID[
      go_map$GO == term
    ]
  )

  if (length(genes_for_term) >= 10) {

    test_go <- term
    test_genes <- genes_for_term

    break
  }
}


if (is.null(test_go)) {
  stop(
    "Could not find a GO term annotated to at least 10 genes."
  )
}


# Keep the technical test reasonably small.

test_genes <- head(
  sort(unique(test_genes)),
  30
)


message(
  "\nTechnical test GO term: ",
  test_go
)

message(
  "Genes selected for smoke test: ",
  length(test_genes)
)


# ---- Determine ontology ------------------------------------------------------

test_ontology <- unique(
  go_map$ONTOLOGY[
    go_map$GO == test_go
  ]
)

test_ontology <- test_ontology[
  !is.na(test_ontology)
]


if (length(test_ontology) == 0) {
  stop(
    "Could not determine ontology for test GO term: ",
    test_go
  )
}


test_ontology <- test_ontology[[1]]


message(
  "Ontology: ",
  test_ontology
)


# ---- Run enrichGO ------------------------------------------------------------

message("\nRunning clusterProfiler::enrichGO()...")


ego <- clusterProfiler::enrichGO(

  gene = test_genes,

  universe = all_gids,

  OrgDb = OrgDb,

  keyType = "GID",

  ont = test_ontology,

  pAdjustMethod = "BH",

  pvalueCutoff = 1,

  qvalueCutoff = 1,

  readable = FALSE
)


# ---- Convert result ----------------------------------------------------------

ego_df <- as.data.frame(ego)


if (nrow(ego_df) == 0) {
  stop(
    "enrichGO() completed but returned no GO terms."
  )
}


# ---- Check whether expected term was recovered ------------------------------

expected_term_found <- test_go %in% ego_df$ID


if (!expected_term_found) {

  warning(
    "Smoke test completed, but the GO term used to construct ",
    "the test gene set was not recovered."
  )

} else {

  message(
    "Expected GO term successfully recovered: ",
    test_go
  )
}


# ---- Save result -------------------------------------------------------------

smoke_output <- file.path(
  dir_output,
  "clusterprofiler_smoke_test.tsv"
)


write.table(
  ego_df,
  smoke_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ---- Show top results --------------------------------------------------------

message("\nTop enrichment results:")

print(
  head(
    ego_df[
      order(ego_df$p.adjust),
    ],
    10
  )
)


# ---- Final validation --------------------------------------------------------

message(
  "\n============================================================\n",
  "CLUSTERPROFILER SMOKE TEST COMPLETE\n",
  "============================================================\n",
  "OrgDb: ",
  pkg_name,
  "\n",
  "Test GO term: ",
  test_go,
  "\n",
  "Ontology: ",
  test_ontology,
  "\n",
  "Test genes: ",
  length(test_genes),
  "\n",
  "Enrichment terms returned: ",
  nrow(ego_df),
  "\n",
  "Expected test term recovered: ",
  expected_term_found,
  "\n",
  "Results: ",
  smoke_output,
  "\n",
  "============================================================"
)

# 02_qc_checks.R
#
# Quality-control checks for the annotation tables prepared by
# 01_prepare_annotation_data.R.
#
# Checks:
#   - gene and GO annotation counts
#   - gene ID overlap between genome and GO annotations
#   - GO annotation coverage
#   - malformed GO identifiers
#   - duplicate gene-GO mappings
#   - GO evidence-code distribution
#
# Outputs:
#   - output/qc_report.txt
#   - output/missing_ids.txt

#Packages#

library(dplyr)
source(here::here("scripts", "00_config.R"))

#Check required inputs#
if (!file.exists(path_gene_info_rds)) {
  stop(
    "gene_info file not found: ",
    path_gene_info_rds,
    "\nRun scripts/01_prepare_annotation_data.R first."
  )
}

if (!file.exists(path_go_final_rds)) {
  stop(
    "go_final file not found: ",
    path_go_final_rds,
    "\nRun scripts/01_prepare_annotation_data.R first."
  )
}

if (!file.exists(path_go_clean_tsv)) {
  stop(
    "Clean GO table not found: ",
    path_go_clean_tsv,
    "\nRun scripts/01_prepare_annotation_data.R first."
  )
}

#Load prepared annotation tables#
gene_info <- readRDS(path_gene_info_rds)
go_final  <- readRDS(path_go_final_rds)

# Read the pre-matching GO table so that unmatched GO identifiers can
# still be detected. go_final itself contains only matched identifiers.
go_all <- read.delim(
  path_go_clean_tsv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

#Basic annotation statistics#
n_genes <- n_distinct(gene_info$GID)

n_go_associations <- nrow(go_all)

n_go_genes <- n_distinct(go_all$GID)

n_go_terms <- n_distinct(go_all$GO)

n_final_associations <- nrow(go_final)

n_final_genes <- n_distinct(go_final$GID)

#Gene ID overlap#
gff_ids <- unique(gene_info$GID)
go_ids  <- unique(go_all$GID)

matched_ids <- intersect(go_ids, gff_ids)
missing_ids <- setdiff(go_ids, gff_ids)

n_matched <- length(matched_ids)
n_missing <- length(missing_ids)

#GO annotation coverage#
annotation_coverage <- if (n_genes > 0) {
  100 * n_final_genes / n_genes
} else {
  NA_real_
}

#Malformed GO IDs
malformed_go <- go_all %>%
  filter(
    is.na(GO) |
      !grepl("^GO:[0-9]{7}$", GO)
  )

n_malformed_go <- nrow(malformed_go)

#Duplicate IDs#
duplicate_go <- go_all %>%
  count(GID, GO, EVIDENCE, name = "n") %>%
  filter(n > 1)

n_duplicate_go <- nrow(duplicate_go)

#Evidence code distribution#
evidence_counts <- go_all %>%
  count(EVIDENCE, name = "n") %>%
  arrange(desc(n))

#Build QC report#
report <- c(
  "============================================================",
  "Coxiella burnetii OrgDb annotation QC",
  "============================================================",
  "",
  "GENOME ANNOTATION",
  sprintf("Unique genes in gene_info              : %d", n_genes),
  "",
  "GO ANNOTATION",
  sprintf("GO associations in source table        : %d", n_go_associations),
  sprintf("Unique GO-annotated gene IDs           : %d", n_go_genes),
  sprintf("Unique GO terms                        : %d", n_go_terms),
  "",
  "GENE ID MATCHING",
  sprintf("GO gene IDs matching gene_info         : %d", n_matched),
  sprintf("GO gene IDs absent from gene_info      : %d", n_missing),
  sprintf("GO associations retained               : %d", n_final_associations),
  "",
  "ANNOTATION COVERAGE",
  sprintf(
    "Genes with >=1 retained GO annotation   : %d",
    n_final_genes
  ),
  sprintf(
    "Genome genes with GO annotation         : %.2f%%",
    annotation_coverage
  ),
  "",
  "GO TABLE QC",
  sprintf("Malformed GO identifiers               : %d", n_malformed_go),
  sprintf("Duplicate GID-GO-EVIDENCE mappings     : %d", n_duplicate_go),
  "",
  "GO EVIDENCE CODES",
  paste0(
    evidence_counts$EVIDENCE,
    " : ",
    evidence_counts$n
  ),
  "",
  "============================================================"
)

#Save QC output#
writeLines(
  report,
  path_qc_report_txt
)

writeLines(
  missing_ids,
  path_missing_ids_txt
)

#print report to console#
message(
  paste(report, collapse = "\n")
)

#Warnings#
if (n_missing > 0) {
  warning(
    n_missing,
    " GO-annotated gene IDs were not found in gene_info. ",
    "See: ",
    path_missing_ids_txt
  )
}

if (n_malformed_go > 0) {
  warning(
    n_malformed_go,
    " malformed GO identifiers were detected."
  )
}

if (n_duplicate_go > 0) {
  warning(
    n_duplicate_go,
    " duplicate GID-GO-EVIDENCE mappings were detected."
  )
}

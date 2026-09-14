
##00_config.R

#Central configuration for the Coxiella burnetii OrgDb build workflow.

#This file defines:
# -project directories
# -input file paths
# -output file paths
# -OrgDb package metadata

# Package required for portable project-relative paths

library(here)

#Project directories

dir_data <- here("data")
dir_output <- here("output")

dir.create(
	dir_output,
	showWarnings = FALSE,
	recursive = TRUE
)

#Input files
#Genome annotation
path_gff <- file.path(
	dir_data,
	"GCF_002094935.1_ASM209493v1_genomic.gff"
)

# Gene Ontology annotation
path_go_raw <- file.path(
  dir_data,
  "GCF_002094935.1_ASM209493v1_gene_ontology.gaf"
)

#Intermediate and output files
path_gene_info_rds <- file.path(
	dir_output,
	"gene_info.rds"
)

path_go_final_rds <- file.path(
  dir_output,
  "go_final.rds"
)

path_go_clean_tsv <- file.path(
  dir_output,
  "coxiella_go_clean.tsv"
)

path_qc_report_txt <- file.path(
  dir_output,
  "qc_report.txt"
)

path_missing_ids_txt <- file.path(
  dir_output,
  "missing_ids.txt"
)

# AnnotationForge will write the generated OrgDb package here.
dir_orgdb_build <- dir_output


# OrgDb package metadata

orgdb_version <- "0.1.0"

orgdb_maintainer <- "Elif Oezer <elif.oezer@uni-wuerzburg.de>"

orgdb_author <- "Elif Oezer"

orgdb_tax_id <- "777"

orgdb_genus <- "Coxiella"

orgdb_species <- "burnetii"

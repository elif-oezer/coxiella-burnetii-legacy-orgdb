#org.Cburnetii.eg.db build pipeline

This repository contains a reproducible workflow for building a Bioconductor-style `OrgDb` annotation package for *Coxiella burnetii* using the **2017 NCBI genome annotation (annotation date: 13 April 2017) associated with assembly `GCA_002094935.1 / ASM209493v1`** and its corresponding Gene Ontology (GO) annotations.

This package is intended for datasets and analyses that use the legacy gene identifier namespace associated with this annotation release. It should not be assumed to be compatible with gene identifiers from newer *C. burnetii* genome annotations.

The resulting package can be used with Bioconductor tools such as `clusterProfiler` for over-representation analysis (ORA) and gene set enrichment analysis (GSEA).


## Compatibility

This OrgDb package is intended for analyses based on the legacy
*Coxiella burnetii* genome annotation and its historical locus-tag
identifier namespace.

The package uses the historical locus tag as the primary gene identifier
(`GID`) because the corresponding RefSeq Gene Ontology annotations are
keyed to this namespace.

It should not be assumed to be compatible with gene identifiers from
newer annotation releases.

A separate OrgDb package for the current genome annotation and current
gene identifiers is planned.

## Repository structure


```
.
├── data/                      # input files and provenance documentation
├── output/                    # generated files (gitignored)
├── scripts/                   # package build workflow
├── .gitignore
└── README.md

````

## Workflow

Genome annotation (GFF)
        +
GO annotation table
        ↓
Prepare gene and GO mappings
        ↓
Quality-control checks
        ↓
Build OrgDb package
        ↓
Install and validate
        ↓
Use in downstream enrichment analysis

## Reproducible workflow

The annotation package is built through a series of reproducible R scripts:

1. `00_setup_environment.R` — installs the required R and Bioconductor dependencies.
2. `01_prepare_annotation_data.R` — imports the genome and GO annotations and prepares the gene and GO mappings.
3. `02_qc_checks.R` — evaluates identifier matching, GO annotation coverage, malformed identifiers, and duplicate mappings.
4. `03_build_orgdb.R` — builds the custom OrgDb package using `AnnotationForge`.
5. `04_install_and_validate_orgdb.R` — installs the generated package and validates its annotation interface.
6. `05_clusterprofiler_smoke_test.R` — performs a functional `clusterProfiler::enrichGO()` smoke test.

Project-specific paths and OrgDb metadata are defined in `00_config.R`.

Generated files and package build artifacts are written to `output/` and are not tracked by Git.

## Input
## Reference annotation

Genome annotation:
- `GCF_002094935.1_ASM209493v1_genomic.gff`

Gene Ontology annotation:
- `GCF_002094935.1_ASM209493v1_gene_ontology.gaf`

Primary OrgDb identifier:
- historical locus tag (`old_locus_tag` in the GFF)

Current locus tags are retained as secondary annotation metadata.

Input provenance and download details are documented in
data/README.md.

## Installation

The OrgDb package is generated locally by the build pipeline and is not
distributed through Bioconductor.

After building the package, it can be installed from the generated package
source directory using:

```r
install.packages(
  "path/to/org.Cburnetii.eg.db",
  repos = NULL,
  type = "source"
)
````
The installed package can then be loaded with:
library(org.Cburnetii.eg.db)

## Basic usage
The package can be queried using AnnotationDbi.
library(AnnotationDbi)
library(org.Cburnetii.eg.db)

keytypes(org.Cburnetii.eg.db)
columns(org.Cburnetii.eg.db)

select(
  org.Cburnetii.eg.db,
  keys = c("BL_example1", "BL_example2"),
  keytype = "GID",
  columns = c("GO", "ONTOLOGY")
)

## GO enrichment with clusterProfiler
The custom OrgDb can be supplied directly to clusterProfiler::enrichGO():

library(clusterProfiler)
library(org.Cburnetii.eg.db)

ego <- enrichGO(
  gene = gene_ids,
  universe = background_gene_ids,
  OrgDb = org.Cburnetii.eg.db,
  keyType = "GID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = FALSE
)





## Status

The current pipeline has been successfully tested from annotation preparation
through OrgDb construction, installation, validation, and a clusterProfiler
GO enrichment smoke test.

The package is currently an initial development release (`0.1.0`).

## License

A license will be added before the first public release.


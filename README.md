#org.Cburnetii.eg.db build pipeline

This repository contains a reproducible workflow for building a 
Bioconductor style "OrgDb" annotation package for *Coxiella burnetii*
using genome annotation and Gene Ontology (GO') data associated with
assembly `GCA_002094935.1 / ASM209493v1`.

The resulting package can be used with Bioconductor tools such as 
`clusterProfiler` for over-representation analysis (ORA) and gene set enrichment analysis (GSEA).


## Repository structure

```text
.
├── data/                      # input files and provenance documentation
├── output/                    # generated files (gitignored)
├── scripts/                   # package build workflow
├── .gitignore
└── README.md


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


## Input

The workflow uses:

genome annotation for assembly GCA_002094935.1 / ASM209493v1
GO annotations derived from
GCF_002094935.1_ASM209493v1_gene_ontology.gaf

Input provenance and download details are documented in
data/README.md.

## Status

This repository is currently being prepared as a reproducible public
release of the OrgDb construction workflow.

## License

A license will be added before the first public release.

```text

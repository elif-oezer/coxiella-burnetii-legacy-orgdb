#org.Cburnetii.eg.db build pipeline

This repository contains a reproducible workflow for building a 
Bioconductor style "OrgDb" annotation package for *Coxiella burnetii*
using genome annotation and Gene Ontology (GO') data associated with
assembly `GCA_002094935.1 / ASM209493v1`.

The resulting package can be used with Bioconductor tools such as 
`clusterProfiler` for over-representation analysis (ORA) and gene set enrichment analysis (GSEA).

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

## Status

The current pipeline has been successfully tested from annotation preparation
through OrgDb construction, installation, validation, and a clusterProfiler
GO enrichment smoke test.

The package is currently an initial development release (`0.1.0`).

## License

A license will be added before the first public release.


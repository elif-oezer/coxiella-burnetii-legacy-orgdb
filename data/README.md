
## Input data

This directory contains the biological annotation files required to build 
the Coxiella burnetii OrgDb package.

The input files themselves are not tracked by Git. Their provenance is
documented here so that the workflow can be reproduced independently.

## Genome annotation

GCA_002094935.1_ASM209493v1_genomic.NZ.gff

Genome annotation file for Coxiella burnetii assembly:

Assembly accession: GCA_002094935.1
Assembly name: ASM209493v1
Source: NCBI
Format: GFF

This file is used to obtain gene identifiers, locus tags, gene symbols, and
gene-product annotations required for construction of the OrgDb gene table.

The file used in the original workflow contains the NZ sequence records from
the assembly annotation.

## Gene Ontology annotation

GCF_002094935.1_ASM209493v1_gene_ontology.gaf

Gene Ontology annotations corresponding to the ASM209493v1 annotation.

Assembly/reference accession: GCF_002094935.1
Assembly name: ASM209493v1
Source: NCBI annotation resources
Format: Gene Association Format (GAF)

The GO annotations used in the OrgDb construction workflow were derived from
this file.

The corresponding mapping used by AnnotationForge::makeOrgPackage() contains
the following three fields:

GID    GO    EVIDENCE

where:

GID is the organism gene/locus identifier;
GO is the Gene Ontology identifier, for example GO:0005524;
EVIDENCE is the GO evidence code.

In the annotation set used for this project, the GO mappings carry the IEA
(Inferred from Electronic Annotation) evidence code.

## Data handling

Raw annotation files are intentionally excluded from Git version control.

The repository .gitignore contains:
data/*
!data/README.md

This keeps the provenance documentation under version control while preventing
downloaded biological annotation files from being committed to the repository.

Generated intermediate files and the resulting OrgDb package are written to
output/, which is also excluded from Git. 

## Expected directory contents

Before running the build workflow, the directory should contain:

data/
├── README.md
├── GCA_002094935.1_ASM209493v1_genomic.NZ.gff
└── GCF_002094935.1_ASM209493v1_gene_ontology.gaf

The build scripts will read these source annotation files and generate the
tables required for construction of the OrgDb package.

# MSH6 Evolution

The `MSH6 Evolution` repository contains code dedicated to the analysis and annotation of the MSH6 protein and its associated domains across the tree of life. This repository provides a  workflow to identify MSH6 orthologs, predict potential reader domains, and annotate these proteins based on various bioinformatics tools and databases.

## Table of Contents
- [Overview](#overview)
- [Purpose](#purpose)
- [Scripts and Functions](#scripts-and-functions)
- [Results Directory Structure](#results-directory-structure)

## Overview

- **Environment Setup**: Utilizes a Conda environment to ensure reproducibility and consistency across different setups.
  
  ```bash
  conda env create -f env.yml
  conda activate tol_reader_repair
  ```

- **InterProScan Installation**: Provides a step-by-step guide to install and set up InterProScan, a tool essential for protein domain annotation.
  
  [Installation Guide](https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html)

- **Scripts & Functions**: Contains a collection of scripts and functions, housed under `code/`, that for tasks such as:
  - Building the necessary directory structure.
  - Downloading and preparing reference genomes.
  - Identifying MSH6 orthologs.
  - Annotating proteins with InterProScan.
  - Parsing results from various bioinformatics tools.

- **Results Directory Structure**: Describes the organization and content of the `results/` directory, which houses the output from various analyses. (Refer to the Results Directory Structure section below)

## Purpose

The primary objective of this repository is to discover and annotate histone reader domains MSH6 orthologs and other repair proteins across different species. The workflow integrates BLAST searches, domain prediction, and annotation of MSH6, other mismatch repair proteins, and their domains. Downstream analyses include measuring sequence conservation of amino acid in histone readers, and phloygenetic contrasts between species with/without histone reader MSH6 fusions.

## Notable `data` Files:

- **`MSH6_arabidopsis_rev_blastp_all.txt`**: Reverse BLASTp results for Arabidopsis MSH6.
- **`MSH6_human_rev_blastp_all.txt`**: Reverse BLASTp results for human MSH6.
- **`tblastn_msh6_blast_dist.txt`**: Calculated distances between MSH6 and domains using tblastn.
- **`msh6_organism_table.txt`**: Summary by organism of predicted MSH6 orthologs annotated for domain state, filtered for ambiguous domain calls.
- **`potential_orthologs_mismatch.txt`**: Predicted non-MSH6 MSH5 proteins.
- **`predicted_msh6_blast.txt`**: Predicted MSH6 orthologs (ortholog to human and A. thaliana MSH6).
- **`predicted_msh6_blast_annotated.txt`**: Predicted MSH6 orthologs annotated with evidence of domain presence/absence and annotation errors.
- **`pwwp_rev_blastp.csv`**: Top hits from blasting PWWP containing proteins against the human genome for each organism.
- **`tudor_rev_blastp.csv`**: Top hits from blasting Tudor containing proteins against the A. thaliana genome for each organism.

# `setup_run.sh` Directory Structure

`setup_run.sh` makes directory structure and runs MSH6 orthology search. Creates `results/` which contains the outputs of various BLAST and annotation processes related to the MSH6 protein and its domains. 

## Subdirectories of `results/` created

### `blastp`
Contains BLASTp results organized by NCBI accession. Each accession directory has subdirectories for different BLASTp queries (e.g., `arabidopsis_msh6`). Inside these, you'll find:
  - GFF files for protein hits.
  - FASTA files for protein hits.
  - Reciprocal BLASTp results against A. thaliana reference.
  - BLASTp query results.
  - Names of proteins in BLASTp query results.

### `tblastn`
Contains tblastn results organized by NCBI accession. Each accession directory has results for different amino acid queries.

### `predicted_msh6`
Contains FASTA files from `human_msh6` and `arabidopsis_msh6` BLASTp that match MSH (including MSH6, MSH1, MSH2, MSH3, etc.). This directory also includes:
  - A combined FASTA file (`all.faa`).
  - Interproscan results in various formats (`gff3`, `json`, `tsv`, `xml`).

## details and structure of results files
```
results/
├── MSH6_arabidopsis_rev_blastp_all.txt 
├── MSH6_human_rev_blastp_all.txt
├── tblastn_msh6_blast_dist.txt # cacluated distance between msh6 and domains with tblastn
├── msh6_organism_table.txt # sumamry by organism: predicted MSH6 orthologs annotated for domain state and filtered for ambiguous domain calls
├── potential_orthologs_mismatch.txt # predicted non-MSH6 MSH5 proteins
├── predicted_msh6_blast.txt # predicted MSH6 orthologs (ortholog to human and A. thaliana MSH6)
├── predicted_msh6_blast_annotated.txt # predicted MSH6 orthologs annotated with evidence of domain presence/absence and annotation error
├── pwwp_rev_blastp.csv #top hits from blasting PWWP containing proteins against human genome for each organism
├── tudor_rev_blastp.csv #top hits from blasting Tudor containing proteins against A. thaliana genome for each organism
├── predicted_msh6 #diretory, see below
├── blastp #diretory, see below
└── tblastn #diretory, see below

results/blastp
├── GCA_000001215.4 #ncbi accession
    ├── arabidopsis_msh6 #blastp query
        ├── gff
        │   ├── AAF49656.2.gff #GFF for protein hit
        │   ├── AAF53392.2.gff
        │   ...
        ├── proteins
        │   ├── AAF49656.2.fasta #fasta for protein hit
        │   ├── AAF53392.2.fasta
        │   ...
        ├── rev_blastp
        │   ├── AAF49656.2.txt #reciprocal blastp against A. thaliana reference format 6
        │   ├── AAF53392.2.txt
        │  ...
        ├── blastp_results.txt #blastp query results format 6
        └── top_hits_ids.txt # names of proteins in blastp query results
    ├── arabidopsis_msh6_tudor
    ├── human_msh6
    ├── human_msh6_pwwp
    └──...
├── GCA_000001735.2
...

results/tblastn
├── GCA_000001215.4
    ├── arabidopsis_msh6 #amino acid query name
    │   └── tblastn_results.txt #tblastn results against genome in format 6
    ├── arabidopsis_msh6_tudor
    ├── human_msh6
    ├── human_msh6_pwwp
    └──...
├── GCA_000001735.2
...

results/predicted_msh6 #all fasta from human_msh6 and arabidopsis_msh6 blastp matching MSH (MSH6 + MSH1, MSH2, MSH3, etc)
├── AAF06041.1.fasta
├── AAF49656.2.fasta
...
├── all.faa # combined fasta
├── all.faa.gff3 # interproscan results
├── all.faa.json
├── all.faa.tsv
└── all.faa.xml

```

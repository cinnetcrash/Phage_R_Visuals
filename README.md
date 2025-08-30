# Phage Genome Visualization in R

This repository contains R scripts to generate circular and linear genome maps of phage genomes based on annotated GFF3 files. The workflow was developed for visualization of *Salmonella* phage genomes annotated with **Bakta** and is suitable for use in comparative genomics or figure preparation for publications.

## Features
- Input: GFF3 annotation files without FASTA sequences (`*.nofasta.gff3`).
- Parsing and normalization of annotation fields.
- Grouping of coding sequences into functional categories:
  - Structural proteins
  - DNA packaging
  - Replication
  - Recombination/Integration
  - Lysis
  - Hypothetical proteins
  - Other
- Visualization:
  - **Linear maps** using [`genoPlotR`](https://cran.r-project.org/package=genoPlotR).
  - **Circular maps** using [`circlize`](https://cran.r-project.org/package=circlize).
  - Separate rings for positive and negative strands.
  - Color-coded functional categories.
  - Genome position grid with tick marks (2 kb intervals, extended to full genome length).
- Optional GC content/skew track if FASTA sequences are available.

## Requirements
- R ≥ 4.1
- Packages: `ape`, `genoPlotR`, `circlize`, (`seqinr` optional for GC skew)

To install missing packages:
```r
install.packages(c("ape","genoPlotR","circlize","seqinr"), dependencies = TRUE)

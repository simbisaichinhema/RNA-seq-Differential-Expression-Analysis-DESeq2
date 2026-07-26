# Chapter 2: Dataset

## GSE201530 — Overview

| Attribute | Details |
|-----------|---------|
| Platform | RNA-seq (bulk) |
| Organism | *Homo sapiens* |
| Tissue | Peripheral blood mononuclear cells (PBMCs) |
| Pathogen | SARS-CoV-2 Omicron variant |
| Sample count | 25+ samples |
| Groups | Healthy control, Vaccinated, Unvaccinated infected |

## Data Source

- Raw count matrix: `data/raw/GSE201530_raw_counts_GRCh38.p13_NCBI.tsv`
- GEO series matrix: `data/raw/GSE201530_series_matrix.txt`
- Source: [NCBI GEO — GSE201530](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201530)

## Comparison Design

This analysis compares two groups:

- **Healthy control** — no infection, no vaccination
- **No vaccination, No prior infection** — unvaccinated individuals acutely infected with Omicron

The comparison reveals the transcriptional signature of acute SARS-CoV-2 Omicron infection in the absence of vaccine-primed immunity.

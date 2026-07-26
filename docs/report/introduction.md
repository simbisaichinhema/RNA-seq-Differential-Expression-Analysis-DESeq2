# Introduction

## Background

RNA sequencing (RNA-seq) has become the gold standard for transcriptomic profiling, enabling researchers to capture a comprehensive snapshot of gene expression across thousands of genes simultaneously. The COVID-19 pandemic has driven an unprecedented wave of transcriptomic studies aimed at understanding the host immune response to SARS-CoV-2 infection and how prior vaccination modulates this response.

GEO dataset **GSE201530** (published in 2025) contains bulk RNA-seq data from peripheral blood mononuclear cells (PBMCs) of individuals infected with the Omicron variant of SARS-CoV-2. The study includes samples from vaccinated and unvaccinated individuals, as well as healthy controls, enabling the investigation of how vaccination history shapes the transcriptomic response to infection.

Prior vaccination has been shown to enhance innate immune priming, leading to faster and more robust antiviral responses upon pathogen encounter. Transcriptomic approaches can reveal the molecular signatures of this trained immunity, identifying specific pathways and gene networks that are differentially activated depending on vaccination status.

## Study Rationale

Understanding which genes and pathways are differentially expressed between vaccinated and unvaccinated individuals upon SARS-CoV-2 infection has several important implications:

- **Vaccine efficacy assessment**: Identifying transcriptional signatures associated with vaccine-mediated protection.
- **Biological insight**: Revealing the molecular mechanisms through which vaccination primes the innate immune system.
- **Biomarker discovery**: Pinpointing genes that could serve as correlates of vaccine-induced protection.
- **Methodological demonstration**: Providing a complete, reproducible bioinformatics workflow for differential expression analysis.

The objective of this analysis is **not** to prove a specific hypothesis but to systematically identify differentially expressed genes and characterize the biological processes and pathways enriched in the dataset.

## Dataset

| Attribute | Value |
|-----------|-------|
| GEO Accession | GSE201530 |
| Organism | *Homo sapiens* |
| Platform | RNA-seq |
| Tissue | Whole blood / PBMC |
| Pathogen | SARS-CoV-2 (Omicron variant) |
| Comparison | Healthy control vs. Unvaccinated infected |

The analysis focuses on the comparison between healthy control samples and individuals with no prior vaccination and no prior infection who were infected with Omicron. This comparison reveals the transcriptional signature of acute unvaccinated SARS-CoV-2 infection in contrast to baseline healthy immune states.

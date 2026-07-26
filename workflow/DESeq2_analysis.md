# DESeq2 Differential Expression Analysis Workflow

## Overview

This document describes the complete workflow for performing differential gene expression analysis on GEO dataset GSE201530 using DESeq2.

## Workflow Steps

### 1. Download GEO Data

**Script**: `scripts/01_download_data.R`

- Accesses GEO via the GEOquery R package
- Downloads the raw count matrix and series matrix
- Saves to `data/raw/`

### 2. Prepare Metadata

**Script**: `scripts/02_prepare_metadata.R`

- Parses the GEO series matrix for sample characteristics
- Filters to the comparison groups of interest:
  - Healthy control
  - No vaccination, No prior infection
- Saves structured metadata to `data/metadata/Metadata.csv`

### 3. DESeq2 Analysis

**Script**: `scripts/03_deseq2_analysis.R`

This is the core analysis script. The workflow within this script is:

1. Load count matrix and metadata
2. Match metadata to count columns
3. Filter to comparison groups
4. Filter lowly expressed genes (≥10 counts in ≥2 samples)
5. Create DESeq2 dataset with design `~ Group`
6. Run DESeq2 (normalization, dispersion estimation, GLM fitting, Wald test)
7. Extract results and map Entrez IDs to gene symbols
8. Save DESeq2 results to `results/tables/DESeq2_Results.csv`
9. Extract significant genes (padj < 0.05, |log2FC| > 1)
10. Export normalized counts
11. Generate PCA, sample distance heatmap, MA plot, volcano plot, and DE gene heatmap

### 4. Visualization

**Script**: `scripts/04_visualization.R`

- Enhanced volcano plot using EnhancedVolcano
- Top 30 DE genes heatmap
- Gene symbol mapping table

### 5. Functional Enrichment

**Script**: `scripts/05_functional_enrichment.R`

- GO enrichment for BP, MF, and CC
- KEGG pathway enrichment
- Dotplots, barplots, enrichment maps, cnetplots, treeplots
- Pathview visualization of top enriched KEGG pathway

## Key Parameters

| Parameter | Value |
|-----------|-------|
| Significance threshold (padj) | 0.05 |
| Fold change threshold (log2FC) | 1.0 |
| Gene expression filter | ≥10 counts in ≥2 samples |
| Multiple testing correction | Benjamini-Hochberg |

## Output Files

| File | Description |
|------|-------------|
| `results/csv/DESeq2_Results.csv` | Full DE results for all genes |
| `results/csv/Significant_Genes.csv` | Significant DE genes |
| `results/csv/Upregulated_Genes.csv` | Upregulated genes only |
| `results/csv/Downregulated_Genes.csv` | Downregulated genes only |
| `results/csv/GO_BP.csv` | GO Biological Process enrichment |
| `results/csv/GO_MF.csv` | GO Molecular Function enrichment |
| `results/csv/GO_CC.csv` | GO Cellular Component enrichment |
| `results/csv/KEGG.csv` | KEGG pathway enrichment |
| `results/csv/Analysis_Summary.csv` | Summary metrics |
| `results/plots/PCA_plot.png/pdf` | PCA plot |
| `results/plots/Sample_Distance_Heatmap.png/pdf` | Sample distance heatmap |
| `results/plots/Heatmap.png/pdf` | Top 50 DE genes heatmap |
| `results/plots/MA_plot.png/pdf` | MA plot |
| `results/plots/Volcano_plot.png/pdf` | Volcano plot |
| `results/plots/GO_BP_Dotplot.png/pdf` | GO BP dotplot |
| `results/plots/GO_BP_Barplot.png/pdf` | GO BP barplot |
| `results/plots/GO_MF_Dotplot.png/pdf` | GO MF dotplot |
| `results/plots/GO_MF_Barplot.png/pdf` | GO MF barplot |
| `results/plots/GO_CC_Dotplot.png/pdf` | GO CC dotplot |
| `results/plots/GO_CC_Barplot.png/pdf` | GO CC barplot |
| `results/plots/KEGG_Dotplot.png/pdf` | KEGG dotplot |
| `results/plots/KEGG_Barplot.png/pdf` | KEGG barplot |
| `results/plots/GO_BP_EnrichmentMap.png/pdf` | GO BP enrichment network |
| `results/plots/GO_BP_Cnetplot.png/pdf` | GO BP gene-concept network |
| `results/plots/GO_BP_Treeplot.png/pdf` | GO BP treeplot |
| `results/plots/KEGG_EnrichmentMap.png/pdf` | KEGG enrichment network |
| `results/plots/KEGG_Cnetplot.png/pdf` | KEGG gene-concept network |
| `results/pathway/hsa05164.png` | KEGG Coronavirus pathway (pathview) |
| `results/pathway/hsa05164.RNAseq.png` | Pathway with DE genes highlighted |

## Reproducibility

All scripts are fully commented and can be re-run from the repository root. Session information is captured in `results/sessionInfo.txt`.

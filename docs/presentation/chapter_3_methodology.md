# Chapter 3: Methodology

## Analysis Pipeline

```
GEO Download → QC → DESeq2 → Visualization → Enrichment → Interpretation
```

## Step-by-Step

1. **Data Download** — Retrieve raw count matrix and metadata from GEO
2. **Quality Control** — Filter lowly expressed genes; inspect PCA and sample distances
3. **DESeq2 Analysis** — Normalize, estimate dispersion, perform Wald test
4. **Visualization** — PCA, MA plot, volcano plot, heatmap
5. **Functional Enrichment** — GO (BP, MF, CC) and KEGG enrichment via clusterProfiler
6. **Pathway Visualization** — pathview for KEGG diagrams
7. **Interpretation** — Biological context and limitations

## Tools Used

| Tool | Purpose |
|------|---------|
| DESeq2 | Differential expression |
| clusterProfiler | GO & KEGG enrichment |
| pathview | KEGG pathway diagrams |
| ggplot2 | Publication figures |
| pheatmap | Sample distances & gene heatmaps |
| EnhancedVolcano | Enhanced volcano plots |

## Key Parameters

- Significance threshold: adjusted p-value < 0.05
- Fold change threshold: |log2FC| > 1
- Gene filter: ≥10 counts in ≥2 samples
- Multiple testing correction: Benjamini-Hochberg

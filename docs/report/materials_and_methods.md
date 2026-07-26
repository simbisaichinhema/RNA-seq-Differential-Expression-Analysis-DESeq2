# Materials and Methods

## Data Source

Raw count data were obtained from the Gene Expression Omnibus (GEO) under accession **GSE201530**. The count matrix (`GSE201530_raw_counts_GRCh38.p13_NCBI.tsv`) was downloaded directly from the GEO repository. Sample metadata was extracted from the GEO series matrix file (`GSE201530_series_matrix.txt`).

## Sample Selection

Samples were filtered into two comparison groups:

1. **Healthy control** — individuals with no COVID-19 infection, no vaccination, and no prior infection.
2. **No vaccination, No prior infection** — individuals infected with the Omicron variant who had no prior vaccination history and no prior SARS-CoV-2 infection.

The comparison aims to capture the intrinsic transcriptional response to acute Omicron infection in the absence of vaccine-primed immunity.

## Quality Control

Quality control was performed using the following metrics:

- **Filtering**: Genes with fewer than 10 counts in fewer than 2 samples were removed prior to analysis. This step removes lowly expressed genes that contribute noise without biological signal.
- **Sample distance**: Pairwise Euclidean distances between variance-stabilized sample expression profiles were computed and visualized as a heatmap to confirm that biological replicates cluster together and that outliers are identifiable.
- **Principal Component Analysis (PCA)**: PCA was performed on the variance-stabilized transformation (VST) of the count data to visualize global sample clustering and confirm separation between groups.

## Differential Expression Analysis

Differential gene expression was performed using **DESeq2** (version 1.40+), a Bioconductor package that models count data using a negative binomial distribution with shrinkage estimation of dispersion and log2 fold changes.

The analysis pipeline is as follows:

1. **Count matrix preparation**: Raw counts were filtered to retain genes with at least 10 counts in at least 2 samples.
2. **DESeqDataSet creation**: A DESeq2 dataset object was constructed using the formula `~ Group`.
3. **Normalization**: DESeq2 applies median-of-ratios normalization to account for differences in library size.
4. **Dispersion estimation**: Gene-wise dispersion estimates are shrunk toward a fitted trend to stabilize variance.
5. **Model fitting**: A negative binomial generalized linear model is fitted for each gene.
6. **Hypothesis testing**: Wald tests are used to test for differential expression between groups.
7. **Result adjustment**: P-values are adjusted for multiple testing using the Benjamini-Hochberg procedure.

### Criteria for Significance

- Adjusted p-value (*padj*) < 0.05
- Absolute log2 fold change > 1

A gene meeting both criteria is classified as significantly differentially expressed.

## Functional Enrichment Analysis

### Gene Ontology (GO) Enrichment

GO enrichment analysis was performed using the **clusterProfiler** package. The top 271 significantly differentially expressed genes (padj < 0.05, |log2FC| > 1) were tested for enrichment in three GO domains:

- **Biological Process (BP)**: Genes annotated with specific biological functions and processes.
- **Molecular Function (MF)**: Genes encoding proteins with specific biochemical activities.
- **Cellular Component (CC)**: Genes localizing to specific subcellular structures.

Enrichment was assessed using the hypergeometric test with Benjamini-Hochberg multiple testing correction (q-value cutoff = 0.05).

### KEGG Pathway Enrichment

KEGG pathway enrichment was also performed using clusterProfiler with the organism set to *hsa* (Homo sapiens). The same gene list and significance thresholds were used.

### Pathway Visualization

The top enriched KEGG pathway (hsa05164 — Coronavirus disease) was visualized using **pathview**, which maps gene-level expression changes onto the KEGG pathway diagram.

## Software and Versions

| Tool | Version |
|------|---------|
| R | ≥ 4.0 |
| DESeq2 | ≥ 1.40 |
| clusterProfiler | ≥ 4.0 |
| enrichplot | ≥ 1.0 |
| pathview | ≥ 1.38 |
| ggplot2 | ≥ 3.4 |
| pheatmap | ≥ 1.0 |
| EnhancedVolcano | ≥ 1.0 |
| GEOquery | ≥ 2.0 |
| Bioconductor | ≥ 3.18 |

## Reproducibility

Every step of this analysis is documented in the `scripts/` directory with fully commented R code. Session information for each analysis step is saved to `results/sessionInfo.txt`. All scripts can be re-run from the repository root to regenerate all results.

## Limitations

1. **Rigid receptor assumption**: DESeq2 treats genes independently and does not model gene-gene interactions.
2. **Comparative design**: The comparison is between healthy controls and unvaccinated infected individuals; vaccinated infected samples are not included in this analysis.
3. **Single time point**: The analysis does not account for temporal dynamics of the immune response across days post-infection.
4. **Bulk RNA-seq**: The use of bulk (rather than single-cell) RNA-seq averages signals across cell types, potentially masking cell-type-specific responses.
5. **Empirical scoring**: DESeq2 uses an empirical Bayes approach for dispersion shrinkage, which may not capture all sources of biological variability.
6. **No explicit solvent or protein–protein interaction modeling**: The enrichment analysis does not model spatial or physical interactions between enriched gene products.

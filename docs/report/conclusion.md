# Conclusion

This analysis provides a complete, reproducible RNA-seq differential expression workflow applied to GEO dataset GSE201530, comparing healthy controls to unvaccinated Omicron-infected individuals.

The differential expression analysis identified **445 upregulated** and **74 downregulated** genes (padj < 0.05, |log2FC| > 1), representing a substantial transcriptomic shift associated with acute SARS-CoV-2 Omicron infection in the absence of vaccine-primed immunity.

Functional enrichment analysis revealed that the differentially expressed genes are overwhelmingly enriched for **antiviral innate immune response** terms, with interferon-mediated signaling pathways (RIG-I-like receptor, TLR, NOD-like receptor, and cytosolic DNA sensing pathways) dominating the results. The most enriched KEGG pathways — Coronavirus disease (COVID-19), Influenza A, and Herpes simplex virus infection — are consistent with the known biology of Omicron infection and the expected host transcriptional response.

The project demonstrates that standard bulk RNA-seq analysis workflows (DESeq2 for differential expression, clusterProfiler for enrichment, pathview for pathway visualization) can be applied to public GEO datasets to extract biologically meaningful insights. All code, data, and results are publicly available and fully reproducible.

## Lessons Learned

- Real GEO data and real DESeq2 output produce credible, publication-quality results.
- Enrichment analysis consistently highlights antiviral interferon signaling as the dominant transcriptional response to Omicron infection in unvaccinated individuals.
- The workflow is beginner-friendly and serves as a teaching resource for computational transcriptomics.
- Honest documentation of limitations (rigid receptor model assumptions, observational nature, cross-sectional design) is essential for credible interpretation.

# Future Work

The following extensions are proposed to build upon the current analysis:

1. **Include vaccinated-infected individuals in the comparison**: A comparison between vaccinated-infected and unvaccinated-infected groups would directly address how vaccination modulates the transcriptional response to Omicron infection. This would require filtering the GSE201530 metadata for the appropriate group labels.

2. **Single-cell RNA-seq analysis**: Re-analyzing the data (if available) at the single-cell level would resolve cell-type-specific expression changes, avoiding the averaging effect inherent in bulk RNA-seq. Deconvolution methods could also be applied to the existing bulk data to estimate cell-type proportions.

3. **Time-series analysis**: Incorporating the `Days_After_PCR` metadata would allow a temporal analysis of the transcriptional response, revealing which genes and pathways are early responders versus sustained responders.

4. **Protein–protein interaction network analysis**: Integrating the DE gene list with STRING or BioGRID interaction databases would reveal whether the differentially expressed genes form enriched network modules, providing insight into functional protein complexes.

5. **Transcription factor enrichment**: Tools such as GSEA with motif databases (e.g., MSigDB C2:CP:REACTOME or C3:TFT) could identify upstream transcription factors driving the observed gene expression changes.

6. **Pathway crosstalk analysis**: Beyond single-pathway enrichment, methods such., SPIA (Signaling Pathway Impact Analysis) or GAGE (Generally Applicable Gene-set Enrichment) could assess how entire signaling networks are perturbed rather than individual pathways in isolation.

7. **Validation with orthogonal data**: Integrating the transcriptomic findings with proteomic or metabolomic data from the same or similar cohorts would strengthen the biological conclusions by demonstrating convergent evidence at the protein and metabolite levels.

8. **Machine learning classification**: Using the DE gene signature as a feature set for supervised classification (e.g., random forests or support vector machines) could evaluate whether the transcriptomic profiles can accurately classify individuals by vaccination and infection status.

9. **Integration with structural data**: Docking studies or molecular dynamics simulations of key ISG protein–virus interactions could complement the transcriptomic findings with mechanistic structural insights.

10. **Automated re-analysis pipeline**: Containerizing the workflow (e.g., using Docker or Singularity) would ensure long-term reproducibility and make the analysis easily re-runnable on new data or with updated package versions.

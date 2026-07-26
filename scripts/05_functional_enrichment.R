############################################################
# 05_functional_enrichment.R
# Purpose: Perform Gene Ontology and KEGG pathway enrichment
#          analysis on differentially expressed genes using
#          clusterProfiler, enrichplot, and pathview.
#
# Required packages: clusterProfiler, org.Hs.eg.db,
#   AnnotationDbi, enrichplot, DOSE, pathview, ggplot2,
#   dplyr, readr
#
# Inputs:
#   results/csv/DESeq2_Results.csv
#
# Outputs:
#   results/csv/GO_BP.csv
#   results/csv/GO_MF.csv
#   results/csv/GO_CC.csv
#   results/csv/KEGG.csv
#   results/csv/Upregulated_Genes.csv
#   results/csv/Downregulated_Genes.csv
#   results/plots/GO_*_Dotplot.png/pdf
#   results/plots/GO_*_Barplot.png/pdf
#   results/plots/GO_*_EnrichmentMap.png/pdf
#   results/plots/GO_*_Cnetplot.png/pdf
#   results/plots/GO_*_Treeplot.png/pdf
#   results/plots/KEGG_*_Dotplot.png/pdf
#   results/plots/KEGG_*_Barplot.png/pdf
#   results/plots/KEGG_*_EnrichmentMap.png/pdf
#   results/plots/KEGG_*_Cnetplot.png/pdf
#   results/pathway/hsa05164.png
#   results/pathway/hsa05164.RNAseq.png
#   results/pathway/hsa05164.xml
#
# Expected runtime: 5–15 minutes
#
# Example execution:
#   Rscript scripts/05_functional_enrichment.R
#
# Error handling: tryCatch blocks for enrichment functions;
#   informative stop messages for missing input files.
############################################################

cat("========================================\n")
cat(" Step 5: Functional Enrichment Analysis\n")
cat("========================================\n\n")

############################################################
# Load Libraries
############################################################

cran_packages <- c("ggplot2", "dplyr", "readr")
bioc_packages <- c(
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
  "enrichplot", "DOSE", "pathview"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  library(pkg, character.only = TRUE)
}

cat("All required packages loaded.\n\n")

############################################################
# Load Utilities
############################################################

source("scripts/utils.R")

############################################################
# Create Output Directories
############################################################

section_header("Setting up output directories")

make_dir("results/go")
make_dir("results/kegg")
make_dir("results/pathway")

############################################################
# Load DESeq2 Results
############################################################

section_header("Loading DESeq2 results")

res_file <- "results/csv/DESeq2_Results.csv"
check_file(res_file)

res <- read.csv(res_file, stringsAsFactors = FALSE)

required_columns <- c("GeneSymbol", "GeneID", "log2FoldChange", "padj")
missing_columns <- setdiff(required_columns, colnames(res))

if (length(missing_columns) > 0) {
  stop(paste0("Missing columns in DESeq2 results: ", paste(missing_columns, collapse = ", ")))
}

res <- res %>%
  filter(
    !is.na(GeneID),
    !is.na(log2FoldChange),
    !is.na(padj)
  )

cat("Total genes loaded:", nrow(res), "\n")

############################################################
# Filter Significant Genes
############################################################

section_header("Filtering significant genes")

upregulated <- res %>%
  filter(padj < 0.05, log2FoldChange > 1)

downregulated <- res %>%
  filter(padj < 0.05, log2FoldChange < -1)

write.csv(upregulated, "results/csv/Upregulated_Genes.csv", row.names = FALSE)
write.csv(downregulated, "results/csv/Downregulated_Genes.csv", row.names = FALSE)

cat("Upregulated genes  :", nrow(upregulated), "\n")
cat("Downregulated genes:", nrow(downregulated), "\n\n")

############################################################
# Prepare Gene List for Enrichment
############################################################

gene_list <- unique(as.character(upregulated$GeneID))

cat("Gene list prepared:", length(gene_list), "unique Entrez IDs\n\n")

############################################################
# GO Biological Process
############################################################

section_header("GO Biological Process enrichment")

go_bp <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

cat("GO BP terms found:", nrow(as.data.frame(go_bp)), "\n")

############################################################
# GO Molecular Function
############################################################

section_header("GO Molecular Function enrichment")

go_mf <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

cat("GO MF terms found:", nrow(as.data.frame(go_mf)), "\n")

############################################################
# GO Cellular Component
############################################################

section_header("GO Cellular Component enrichment")

go_cc <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

cat("GO CC terms found:", nrow(as.data.frame(go_cc)), "\n")

############################################################
# KEGG Pathway Enrichment
############################################################

section_header("KEGG pathway enrichment")

kegg <- enrichKEGG(
  gene = gene_list,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

cat("KEGG pathways found:", nrow(as.data.frame(kegg)), "\n")

############################################################
# Save Enrichment Tables
############################################################

section_header("Saving enrichment tables")

write.csv(as.data.frame(go_bp), "results/csv/GO_BP.csv", row.names = FALSE)
write.csv(as.data.frame(go_mf), "results/csv/GO_MF.csv", row.names = FALSE)
write.csv(as.data.frame(go_cc), "results/csv/GO_CC.csv", row.names = FALSE)
write.csv(as.data.frame(kegg), "results/csv/KEGG.csv", row.names = FALSE)

cat("Enrichment tables saved to results/csv/\n\n")

############################################################
# Plot Saving Function
############################################################

save_plot <- function(filename, plot_object) {
  png(filename = filename, width = 2400, height = 1800, res = 300)
  print(plot_object)
  dev.off()

  pdf(file = sub("\\.png$", ".pdf", filename), width = 12, height = 9)
  print(plot_object)
  dev.off()
}

############################################################
# GO Dotplots
############################################################

section_header("GO dotplots")

if (nrow(as.data.frame(go_bp)) > 0) {
  save_plot("results/plots/GO_BP_Dotplot.png", dotplot(go_bp, showCategory = 20, title = "GO Biological Process"))
}
if (nrow(as.data.frame(go_mf)) > 0) {
  save_plot("results/plots/GO_MF_Dotplot.png", dotplot(go_mf, showCategory = 20, title = "GO Molecular Function"))
}
if (nrow(as.data.frame(go_cc)) > 0) {
  save_plot("results/plots/GO_CC_Dotplot.png", dotplot(go_cc, showCategory = 20, title = "GO Cellular Component"))
}

cat("GO dotplots saved.\n")

############################################################
# GO Barplots
############################################################

section_header("GO barplots")

if (nrow(as.data.frame(go_bp)) > 0) {
  save_plot("results/plots/GO_BP_Barplot.png", barplot(go_bp, showCategory = 20, title = "GO Biological Process"))
}
if (nrow(as.data.frame(go_mf)) > 0) {
  save_plot("results/plots/GO_MF_Barplot.png", barplot(go_mf, showCategory = 20, title = "GO Molecular Function"))
}
if (nrow(as.data.frame(go_cc)) > 0) {
  save_plot("results/plots/GO_CC_Barplot.png", barplot(go_cc, showCategory = 20, title = "GO Cellular Component"))
}

cat("GO barplots saved.\n")

############################################################
# KEGG Plots
############################################################

section_header("KEGG plots")

if (nrow(as.data.frame(kegg)) > 0) {
  save_plot("results/plots/KEGG_Dotplot.png", dotplot(kegg, showCategory = 20, title = "KEGG Pathway Enrichment"))
  save_plot("results/plots/KEGG_Barplot.png", barplot(kegg, showCategory = 20, title = "KEGG Pathway Enrichment"))
}

cat("KEGG plots saved.\n")

############################################################
# GO Enrichment Map + Network Plots
############################################################

section_header("GO enrichment network plots")

bp_sim <- tryCatch(pairwise_termsim(go_bp), error = function(e) NULL)

if (!is.null(bp_sim)) {
  save_plot("results/plots/GO_BP_EnrichmentMap.png", emapplot(bp_sim, showCategory = 20))
  save_plot("results/plots/GO_BP_Cnetplot.png", cnetplot(bp_sim, showCategory = 10))
  save_plot("results/plots/GO_BP_Treeplot.png", treeplot(bp_sim))
  cat("GO BP enrichment network plots saved.\n")
} else {
  cat("GO BP enrichment map skipped.\n")
}

############################################################
# KEGG Network Plots
############################################################

section_header("KEGG network plots")

kegg_sim <- tryCatch(pairwise_termsim(kegg), error = function(e) NULL)

if (!is.null(kegg_sim)) {
  save_plot("results/plots/KEGG_EnrichmentMap.png", emapplot(kegg_sim, showCategory = 20))
  save_plot("results/plots/KEGG_Cnetplot.png", cnetplot(kegg_sim, showCategory = 10))
  cat("KEGG enrichment network plots saved.\n")
} else {
  cat("KEGG enrichment map skipped.\n")
}

############################################################
# Pathview
############################################################

section_header("Generating KEGG pathway diagram")

fold_change <- upregulated$log2FoldChange
names(fold_change) <- as.character(upregulated$GeneID)

if (nrow(as.data.frame(kegg)) > 0) {
  top_pathway <- as.data.frame(kegg)$ID[1]

  tryCatch(
    {
      pathview(
        gene.data = fold_change,
        pathway.id = top_pathway,
        species = "hsa",
        out.suffix = "RNAseq"
      )
      cat("Pathview diagram generated for:", top_pathway, "\n")
    },
    error = function(e) {
      cat("Pathview failed:", e$message, "\n")
    }
  )
} else {
  cat("No enriched KEGG pathway found — pathview skipped.\n")
}

############################################################
# Save Session Info
############################################################

session_info("results/sessionInfo.txt")

############################################################
# Analysis Summary
############################################################

section_header("Analysis Complete")

summary_report <- data.frame(
  Metric = c(
    "Total Genes",
    "Upregulated Genes",
    "Downregulated Genes",
    "GO Biological Process",
    "GO Molecular Function",
    "GO Cellular Component",
    "KEGG Pathways"
  ),
  Value = c(
    nrow(res),
    nrow(upregulated),
    nrow(downregulated),
    nrow(as.data.frame(go_bp)),
    nrow(as.data.frame(go_mf)),
    nrow(as.data.frame(go_cc)),
    nrow(as.data.frame(kegg))
  )
)

write.csv(summary_report, "results/csv/Analysis_Summary.csv", row.names = FALSE)

cat("\nSummary Report:\n")
print(summary_report, row.names = FALSE)

cat("\nAll functional enrichment results saved.\n")
cat("  results/csv/GO_BP.csv\n")
cat("  results/csv/GO_MF.csv\n")
cat("  results/csv/GO_CC.csv\n")
cat("  results/csv/KEGG.csv\n")
cat("  results/plots/*\n")
cat("  results/pathway/*\n")

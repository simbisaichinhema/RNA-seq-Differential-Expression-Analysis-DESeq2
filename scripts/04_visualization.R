############################################################
# 04_visualization.R
# Purpose: Generate supplementary and publication-quality
#          visualizations from the DESeq2 results.
#
# This script produces figures that complement the plots
# generated in 03_deseq2_analysis.R, including enhanced
# volcano plots and additional heatmaps.
#
# Required packages: DESeq2, ggplot2, pheatmap, RColorBrewer,
#                    org.Hs.eg.db, AnnotationDbi, ggrepel,
#                    EnhancedVolcano
#
# Inputs:
#   results/csv/DESeq2_Results.csv
#   results/csv/Normalized_Counts.csv
#   data/metadata/Metadata.csv
#
# Outputs:
#   results/plots/* (PNG + PDF)
#   results/csv/* (CSV tables)
#
# Expected runtime: 2–5 minutes
#
# Example execution:
#   Rscript scripts/04_visualization.R
#
# Error handling: Checks for required input files and
#   column names; stops with informative messages.
############################################################

cat("========================================\n")
cat(" Step 4: Additional Visualizations\n")
cat("========================================\n\n")

############################################################
# Load Libraries
############################################################

required_packages <- c(
  "DESeq2", "ggplot2", "pheatmap", "RColorBrewer",
  "org.Hs.eg.db", "AnnotationDbi", "ggrepel", "EnhancedVolcano"
)

for (pkg in required_packages) {
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
# Load Data
############################################################

section_header("Loading results")

res_df <- read.csv("results/csv/DESeq2_Results.csv", stringsAsFactors = FALSE)
metadata <- read.csv("data/metadata/Metadata.csv", stringsAsFactors = FALSE)

validate_results(res_df)

# Convert Entrez IDs to gene symbols if not already present
if (!"GeneSymbol" %in% colnames(res_df)) {
  message("Mapping Entrez IDs to gene symbols...")
  res_df$GeneSymbol <- mapIds(
    org.Hs.eg.db,
    keys = res_df$GeneID,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
}

cat("Loaded", nrow(res_df), "genes\n")
cat("Significant genes (padj < 0.05, |log2FC| > 1):",
    sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, na.rm = TRUE), "\n\n")

############################################################
# Enhanced Volcano Plot
############################################################

section_header("Generating enhanced volcano plot")

p_enhanced <- EnhancedVolcano(
  res_df,
  lab = res_df$GeneSymbol,
  x = "log2FoldChange",
  y = "padj",
  title = "GSE201530: Enhanced Volcano Plot",
  subtitle = "Healthy control vs. Unvaccinated Omicron-infected",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 3.0,
  colAlpha = 0.7,
  legendPosition = "right",
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  maxoverlaps = Inf
)

ggsave("results/plots/Enhanced_Volcano.png", p_enhanced, width = 10, height = 8, dpi = 300)
ggsave("results/plots/Enhanced_Volcano.pdf", p_enhanced, width = 10, height = 8)

cat("Enhanced volcano plot saved.\n")

############################################################
# Top 30 DE Genes Heatmap (with gene symbols)
############################################################

section_header("Generating top 30 DE genes heatmap")

sig_genes <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
top30 <- head(sig_genes[order(sig_genes$pvalue), ], 30)

counts_norm <- read.csv("results/csv/Normalized_Counts.csv", row.names = 1, check.names = FALSE)

# Subset to top 30 genes
mat_top30 <- counts_norm[rownames(counts_norm) %in% top30$GeneID, ]

# Use gene symbols as row names
rownames(mat_top30) <- top30$GeneSymbol[match(rownames(mat_top30), top30$GeneID)]

# Z-score normalization per gene
mat_top30 <- t(scale(t(mat_top30)))

annotation_col <- data.frame(Group = metadata$Group)
rownames(annotation_col) <- metadata$Sample

pheatmap(
  mat_top30,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 7,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  filename = "results/plots/Top30_DE_Heatmap.png",
  width = 10,
  height = 12
)
ggsave("results/plots/Top30_DE_Heatmap.pdf",
       pheatmap(mat_top30, scale = "row", annotation_col = annotation_col, show_rownames = TRUE, fontsize_row = 7, cluster_rows = TRUE, cluster_cols = TRUE, color = colorRampPalette(c("navy", "white", "firebrick3"))(100)),
       width = 10, height = 12)

cat("Top 30 DE genes heatmap saved.\n\n")

############################################################
# Gene ID to Symbol mapping helper table
############################################################

section_header("Saving gene mapping table")

gene_map <- data.frame(
  GeneID = res_df$GeneID,
  GeneSymbol = res_df$GeneSymbol,
  stringsAsFactors = FALSE
)
gene_map <- gene_map[!duplicated(gene_map$GeneID), ]

write.csv(gene_map, "results/csv/Gene_Symbol_Mapping.csv", row.names = FALSE)
cat("Gene mapping table saved.\n\n")

############################################################
# Summary
############################################################

section_header("Step 4 Complete")
cat("Supplementary visualizations saved to results/plots/\n")

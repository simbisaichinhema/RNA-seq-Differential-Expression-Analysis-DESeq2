############################################################
# 03_deseq2_analysis.R
# Purpose: Perform differential gene expression analysis
#          using DESeq2 on GSE201530 RNA-seq data.
#
# Comparison: Healthy control vs. No vaccination,
#             No prior infection (unvaccinated infected).
#
# Required packages: DESeq2, ggplot2, pheatmap, RColorBrewer,
#                    genefilter, org.Hs.eg.db, AnnotationDbi,
#                    ggrepel
#
# Inputs:
#   data/raw/GSE201530_raw_counts_GRCh38.p13_NCBI.tsv
#   data/metadata/Metadata.csv
#
# Outputs:
#   results/tables/DESeq2_Results.csv
#   results/tables/Significant_Genes.csv
#   results/tables/Normalized_Counts.csv
#   results/tables/Volcano_Labeled_Genes.csv
#   results/tables/Metadata.csv (filtered)
#   results/plots/PCA_plot.png / .pdf
#   results/plots/Sample_Distance_Heatmap.png / .pdf
#   results/plots/MA_plot.png / .pdf
#   results/plots/Volcano_plot.png / .pdf
#   results/plots/Heatmap.png / .pdf
#
# Expected runtime: 3–8 minutes
#
# Example execution:
#   Rscript scripts/03_deseq2_analysis.R
#
# Error handling: Stops with informative messages for
#   mismatched metadata, empty results, or failed DESeq2 runs.
############################################################

cat("========================================\n")
cat(" Step 3: DESeq2 Differential Expression\n")
cat("========================================\n\n")

############################################################
# Load Libraries
############################################################

required_packages <- c(
  "DESeq2", "ggplot2", "pheatmap", "RColorBrewer",
  "genefilter", "org.Hs.eg.db", "AnnotationDbi", "ggrepel"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Required package not installed: ", pkg, ". Install with BiocManager::install('", pkg, "')."))
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

make_dir("results")
make_dir("results/plots")
make_dir("results/tables")
make_dir("results/csv")

############################################################
# Load Count Matrix
############################################################

section_header("Loading count matrix")

counts_file <- "data/raw/GSE201530_raw_counts_GRCh38.p13_NCBI.tsv"
check_file(counts_file)

counts <- read.delim(counts_file, check.names = FALSE)
cat("Loaded count matrix:", nrow(counts), "genes x", ncol(counts), "columns\n")

############################################################
# Load Metadata
############################################################

section_header("Loading metadata")

metadata_file <- "data/metadata/Metadata.csv"
check_file(metadata_file)

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
cat("Loaded metadata:", nrow(metadata), "samples\n")

############################################################
# Match Metadata to Count Matrix
############################################################

section_header("Matching metadata to count matrix")

# Exclude the GeneID column from sample names
sample_cols <- setdiff(colnames(counts), "GeneID")

metadata <- metadata[metadata$GSM %in% sample_cols, ]
metadata <- metadata[match(sample_cols, metadata$GSM), ]
rownames(metadata) <- metadata$GSM

if (!all(rownames(metadata) == sample_cols)) {
  stop("Metadata and count matrix sample names do not match.")
}

cat("All", ncol(counts) - 1, "samples matched successfully.\n\n")

############################################################
# Filter Comparison Groups
############################################################

section_header("Filtering comparison groups")

selected_groups <- c("Healthy control", "No vaccination, No prior infection")
metadata_filtered <- metadata[metadata$Group %in% selected_groups, ]
metadata_filtered$Group <- factor(
  metadata_filtered$Group,
  levels = c("Healthy control", "No vaccination, No prior infection")
)

cat("Control group:", sum(metadata_filtered$Group == "Healthy control"), "samples\n")
cat("Unvaccinated infected:", sum(metadata_filtered$Group == "No vaccination, No prior infection"), "samples\n\n")

############################################################
# Filter Count Matrix
############################################################

section_header("Filtering count matrix")

counts_filtered <- counts[, c("GeneID", rownames(metadata_filtered)), drop = FALSE]
count_df <- counts_filtered
counts_matrix <- count_df[, -which(names(count_df) == "GeneID"), drop = FALSE]
rownames(counts_matrix) <- count_df$GeneID
counts_matrix <- as.matrix(counts_matrix)
mode(counts_matrix) <- "integer"

cat("Count matrix dimensions:", nrow(counts_matrix), "genes x", ncol(counts_matrix), "samples\n\n")

############################################################
# Filter Lowly Expressed Genes
############################################################

section_header("Filtering lowly expressed genes")

keep <- rowSums(counts_matrix >= 10) >= 2
counts_filtered <- counts_matrix[keep, ]
cat("Retained", nrow(counts_filtered), "genes after filtering (min 10 counts in >= 2 samples)\n\n")

############################################################
# Create DESeq2 Dataset
############################################################

section_header("Creating DESeq2 dataset")

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata_filtered,
  design = ~ Group
)

cat("DESeq2 dataset created with", nrow(dds), "genes and", ncol(dds), "samples\n")
cat("Design formula: ~ Group\n\n")

############################################################
# Run DESeq2 Differential Expression
############################################################

section_header("Running DESeq2 normalization and dispersion estimation")

dds <- DESeq(dds)

cat("DESeq2 completed.\n")
cat("Dispersion:", summary(dds)$dispersion, "\n\n")

############################################################
# Extract and Format Results
############################################################

section_header("Extracting differential expression results")

res <- results(dds)
res <- res[order(res$padj), ]
res_df <- as.data.frame(res)

############################################################
# Convert Entrez IDs to Gene Symbols
############################################################

cat("Mapping Entrez IDs to gene symbols...\n")

res_df$GeneID <- rownames(res_df)

res_df$GeneSymbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res_df),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

res_df <- res_df[, c("GeneSymbol", "GeneID", setdiff(colnames(res_df), c("GeneSymbol", "GeneID")))]

cat("Done.\n")
cat("Total results:", nrow(res_df), "genes\n")
cat("Genes with adjusted p-value < 0.05:", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("Significant genes (padj < 0.05, |log2FC| > 1):", sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, na.rm = TRUE), "\n\n")

############################################################
# Save DESeq2 Results
############################################################

section_header("Saving results")

write.csv(res_df, "results/tables/DESeq2_Results.csv", row.names = FALSE)
write.csv(res_df, "results/csv/DESeq2_Results.csv", row.names = FALSE)
cat("DESeq2 results saved.\n")

############################################################
# Extract Significant Genes
############################################################

sig_genes <- subset(
  res_df,
  padj < 0.05 & abs(log2FoldChange) > 1
)

write.csv(sig_genes, "results/tables/Significant_Genes.csv", row.names = FALSE)
write.csv(sig_genes, "results/csv/Significant_Genes.csv", row.names = FALSE)
cat("Significant genes saved:", nrow(sig_genes), "genes\n")

############################################################
# Count Up/Downregulated
############################################################

n_up <- sum(sig_genes$log2FoldChange > 1, na.rm = TRUE)
n_down <- sum(sig_genes$log2FoldChange < -1, na.rm = TRUE)

cat("Upregulated:", n_up, "\n")
cat("Downregulated:", n_down, "\n\n")

############################################################
# Export Normalized Counts
############################################################

section_header("Exporting normalized counts")

normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, "results/tables/Normalized_Counts.csv", row.names = TRUE)
write.csv(normalized_counts, "results/csv/Normalized_Counts.csv", row.names = TRUE)
cat("Normalized counts saved.\n")

############################################################
# Variance Stabilizing Transformation
############################################################

section_header("Variance stabilizing transformation")

vsd <- vst(dds, blind = FALSE)
cat("VST completed.\n\n")

############################################################
# PCA Plot
############################################################

section_header("Generating PCA plot")

pcaData <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(
  pcaData,
  aes(PC1, PC2, color = Group)
) +
  geom_point(size = 4) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  ggtitle("PCA Plot — GSE201530") +
  theme_classic(base_size = 14)

ggsave("results/plots/PCA_plot.png", p_pca, width = 8, height = 6, dpi = 300)
ggsave("results/plots/PCA_plot.pdf", p_pca, width = 8, height = 6)
ggsave("results/csv/PCA_plot.png", p_pca, width = 8, height = 6, dpi = 300)
cat("PCA plot saved.\n")

############################################################
# Sample Distance Heatmap
############################################################

section_header("Generating sample distance heatmap")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- metadata_filtered$Sample
colnames(sampleDistMatrix) <- metadata_filtered$Sample

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
  filename = "results/plots/Sample_Distance_Heatmap.png",
  width = 10,
  height = 10
)
ggsave("results/plots/Sample_Distance_Heatmap.pdf", pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)), width = 10, height = 10)

cat("Sample distance heatmap saved.\n")

############################################################
# MA Plot
############################################################

section_header("Generating MA plot")

png("results/plots/MA_plot.png", width = 1800, height = 1500, res = 300)
plotMA(res, ylim = c(-6, 6))
dev.off()
ggsave("results/plots/MA_plot.pdf", plotMA(res, ylim = c(-6, 6)), width = 10, height = 8)

cat("MA plot saved.\n")

############################################################
# Volcano Plot
############################################################

section_header("Generating volcano plot")

TOP_LABELS <- 15
SIG_PADJ <- 0.05
SIG_LFC <- 1

volcano <- res_df
volcano$Gene <- volcano$GeneSymbol
volcano <- volcano[!is.na(volcano$padj), ]

volcano$Significance <- "Not Significant"
volcano$Significance[volcano$padj < SIG_PADJ & volcano$log2FoldChange > SIG_LFC] <- "Upregulated"
volcano$Significance[volcano$padj < SIG_PADJ & volcano$log2FoldChange < -SIG_LFC] <- "Downregulated"

volcano$rank <- order(volcano$padj, decreasing = FALSE)
top_genes <- head(volcano[order(volcano$padj), ], TOP_LABELS)
top_genes <- top_genes[!is.na(top_genes$GeneSymbol) & top_genes$GeneSymbol != "", ]

n_up_v <- sum(volcano$Significance == "Upregulated")
n_down_v <- sum(volcano$Significance == "Downregulated")
n_ns <- sum(volcano$Significance == "Not Significant")

cat("Upregulated:", n_up_v, "| Downregulated:", n_down_v, "| Not significant:", n_ns, "\n")

p_volcano <- ggplot(
  volcano,
  aes(x = log2FoldChange, y = -log10(padj), color = Significance)
) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = c(-SIG_LFC, SIG_LFC), linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = -log10(SIG_PADJ), linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_text_repel(
    data = top_genes,
    aes(label = Gene),
    size = 3, fontface = "bold", max.overlaps = Inf,
    box.padding = 0.4, point.padding = 0.3,
    segment.color = "grey50", segment.size = 0.3,
    arrow = arrow(length = unit(0.01, "npc")),
    show.legend = FALSE
  ) +
  scale_color_manual(
    name = "Significance",
    values = c(
      "Downregulated" = "#2166AC",
      "Not Significant" = "grey70",
      "Upregulated" = "#B2182B"
    )
  ) +
  labs(
    title = "Volcano Plot of Differential Expression",
    subtitle = paste0("n = ", nrow(volcano), " | Up: ", n_up_v, " | Down: ", n_down_v),
    x = expression(log[2] ~ "Fold Change"),
    y = expression(-log[10] ~ "Adjusted P-value")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("results/plots/Volcano_plot.png", p_volcano, width = 8, height = 6, dpi = 300)
ggsave("results/plots/Volcano_plot.pdf", p_volcano, width = 8, height = 6)

write.csv(
  top_genes[, c("GeneSymbol", "GeneID", "log2FoldChange", "pvalue", "padj", "Significance")],
  "results/tables/Volcano_Labeled_Genes.csv",
  row.names = FALSE
)
write.csv(
  top_genes[, c("GeneSymbol", "GeneID", "log2FoldChange", "pvalue", "padj", "Significance")],
  "results/csv/Volcano_Labeled_Genes.csv",
  row.names = FALSE
)

cat("Volcano plot saved.\n\n")

############################################################
# Heatmap of Top 50 Differentially Expressed Genes
############################################################

section_header("Generating DE gene heatmap")

top50 <- head(order(res$pvalue), 50)
mat <- assay(vsd)[top50, ]

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(mat),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

rownames(mat) <- ifelse(is.na(gene_symbols), rownames(mat), gene_symbols)
mat <- mat - rowMeans(mat)

annotation <- data.frame(Group = metadata_filtered$Group)
rownames(annotation) <- rownames(metadata_filtered)

pheatmap(
  mat,
  scale = "row",
  annotation_col = annotation,
  show_rownames = TRUE,
  fontsize_row = 8,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  filename = "results/plots/Heatmap.png"
)
ggsave("results/plots/Heatmap.pdf", pheatmap(mat, scale = "row", annotation_col = annotation, show_rownames = TRUE, fontsize_row = 8, cluster_rows = TRUE, cluster_cols = TRUE, color = colorRampPalette(c("navy", "white", "firebrick3"))(100)), width = 12, height = 10)

cat("Heatmap saved.\n\n")

############################################################
# Save Session Info
############################################################

session_info("results/sessionInfo.txt")

############################################################
# Summary
############################################################

section_header("Analysis Complete")
cat("Results saved to:\n")
cat("  results/tables/\n")
cat("  results/csv/\n")
cat("  results/plots/\n")
cat("\nSummary:\n")
cat("  Total genes analysed:", nrow(res_df), "\n")
cat("  Upregulated:", n_up, "\n")
cat("  Downregulated:", n_down, "\n")
cat("  Significant (padj < 0.05, |log2FC| > 1):", nrow(sig_genes), "\n")

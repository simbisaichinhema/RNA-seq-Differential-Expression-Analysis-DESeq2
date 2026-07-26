############################################################
# RNA-Seq Differential Expression Analysis Pipeline
#
# Dataset : GSE201530
# Organism: Homo sapiens
# Platform: RNA-Seq
#
# Author : Simbisai Chinhema
# GitHub :
#
############################################################

rm(list = ls())

cat("========================================\n")
cat(" GSE201530 RNA-Seq Analysis Pipeline\n")
cat("========================================\n\n")

############################################################
# Load Libraries
############################################################

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(org.Hs.eg.db)
library(AnnotationDbi)
########### V2 UPGRADE ###########
library(ggrepel)
######## END V2 UPGRADE ########

############################################################
# Create Output Directories
############################################################

dir.create("results", showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

############################################################
# Read Count Matrix
############################################################

cat("Reading Count Matrix...\n")

counts <- read.delim("GSE201530_raw_counts_GRCh38.p13_NCBI.tsv", check.names = FALSE)

cat("Done.\n\n")

############################################################
# Read GEO Series Matrix
############################################################

cat("Reading GEO Metadata...\n")

geo <- readLines("GSE201530_series_matrix.txt")

cat("Done.\n\n")

############################################################
# Extract GSM IDs
############################################################

gsm <- sub("!Sample_geo_accession = ", "", geo[grep("^!Sample_geo_accession", geo)])
gsm <- strsplit(gsm, "\t")[[1]]
gsm <- gsub('"', "", gsm)

############################################################
# Extract Sample Titles
############################################################

titles <- sub("!Sample_title = ", "", geo[grep("^!Sample_title", geo)])
sample <- strsplit(titles, "\t")[[1]]
sample <- gsub('"', "", sample)

############################################################
# Parse Metadata
############################################################

char_lines <- geo[grep("^!Sample_characteristics_ch1", geo)]

extract_metadata <- function(line) {
    values <- strsplit(line, "\t")[[1]][-1]
    values <- gsub('"', "", values)
    field <- sub(":.*", "", values[1])
    field <- trimws(field)
    values <- sub("^[^:]+: ", "", values)
    list(field = field, values = values)
}

metadata_list <- lapply(
    char_lines,
    extract_metadata
)

names(metadata_list) <- sapply(
    metadata_list,
    function(x) x$field
)

############################################################
# Build Metadata Table
############################################################

# Pad any missing metadata values with NA so every field has
# exactly length(gsm) entries. Missing characteristics for a
# sample will appear as NA in the final table.
n_samples <- length(gsm)
for (.nm in names(metadata_list)) {
    .n <- length(metadata_list[[.nm]]$values)
    if (.n < n_samples) {
        metadata_list[[.nm]]$values <- c(
            metadata_list[[.nm]]$values,
            rep(NA_character_, n_samples - .n)
        )
    }
}
rm(.nm, .n, n_samples)

metadata <- data.frame(
    GSM = gsm,
    Sample = sample,
    Gender = metadata_list[["gender"]]$values,
    Age = as.numeric(metadata_list[["age"]]$values),
    Group = metadata_list[["group (by covid-19 vaccination, prior infection)"]]$values,
    Omicron_Lineage = metadata_list[["omicron sublineage"]]$values,
    Days_After_PCR = metadata_list[["days after positive pcr results"]]$values,
    Disease_State = metadata_list[["disease state"]]$values,
    Country = metadata_list[["geographical location"]]$values,
    Cell_Type = metadata_list[["cell type"]]$values,
    stringsAsFactors = FALSE
)

############################################################
# Save Metadata
############################################################

write.csv(metadata, "results/tables/Metadata.csv", row.names = FALSE)

############################################################
# Match Metadata to Count Matrix
############################################################

cat("Matching Metadata to Count Matrix...\n")

metadata <- metadata[metadata$GSM %in% colnames(counts), ]
metadata <- metadata[match(colnames(counts)[-1], metadata$GSM), ]
rownames(metadata) <- metadata$GSM

if (!all(rownames(metadata) == colnames(counts)[-1])) {
    stop("Metadata and count matrix do not match.")
}

cat("Done.\n\n")

############################################################
# Select Comparison Groups
############################################################

cat("Selecting Comparison Groups...\n")

selected_groups <- c(
    "Healthy control",
    "No vaccination, No prior infection"
)

metadata_filtered <- metadata[metadata$Group %in% selected_groups, ]

metadata_filtered$Group <- factor(
    metadata_filtered$Group,
    levels = c("Healthy control", "No vaccination, No prior infection")
)

cat("Done.\n\n")

############################################################
# Filter Count Matrix
############################################################

cat("Filtering Count Matrix...\n")

counts_filtered <- counts[, c("GeneID", rownames(metadata_filtered)), drop = FALSE]

cat("Done.\n\n")

############################################################
# Prepare Count Matrix
############################################################

cat("Preparing Count Matrix...\n")

########### V2 UPGRADE ###########
# Preserve the original count table (data.frame) for GeneID access
# Maintain count_df (data.frame) separately from counts (matrix for DESeq2)
count_df <- counts_filtered

# Create the count matrix for DESeq2 (drops GeneID column from the matrix)
counts <- count_df[, -which(names(count_df) == "GeneID"), drop = FALSE]
rownames(counts) <- count_df$GeneID
counts <- as.matrix(counts)
mode(counts) <- "integer"
########### END V2 UPGRADE ##########

############################################################
# Filter Lowly Expressed Genes
############################################################

cat("Filtering Lowly Expressed Genes...\n")

keep <- rowSums(counts >= 10) >= 2
counts <- counts[keep, ]

cat(nrow(counts), "genes retained after filtering.\n\n")

############################################################
# Create DESeq2 Dataset
############################################################

cat("Creating DESeq2 Dataset...\n")

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata_filtered,
    design = ~ Group
)

cat("Done.\n\n")

############################################################
# Differential Expression Analysis
############################################################

cat("Running DESeq2...\n")

dds <- DESeq(dds)

cat("Done.\n\n")

############################################################
# Extract Results
############################################################

cat("Extracting Results...\n")

res <- results(dds)
res <- res[order(res$padj), ]
res_df <- as.data.frame(res)

############################################################
# Convert Entrez IDs to Gene Symbols
############################################################

res_df$GeneID <- rownames(res_df)

res_df$GeneSymbol <- mapIds(
    org.Hs.eg.db,
    keys = rownames(res_df),
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
)

res_df <- res_df[, c("GeneSymbol", "GeneID", setdiff(colnames(res_df), c("GeneSymbol", "GeneID")))]

cat("Done.\n\n")

############################################################
# Save DESeq2 Results
############################################################

write.csv(
    res_df,
    "results/tables/DESeq2_Results.csv"
)

############################################################
# Extract Significant Genes
############################################################

sig_genes <- subset(
    res_df,
    padj < 0.05 &
    abs(log2FoldChange) > 1
)

write.csv(
    sig_genes,
    "results/tables/Significant_Genes.csv"
)

############################################################
# Export Normalized Counts
############################################################

normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, "results/tables/Normalized_Counts.csv")

############################################################
# Variance Stabilizing Transformation
############################################################

cat("Performing Variance Stabilizing Transformation...\n")

vsd <- vst(dds, blind = FALSE)

cat("Done.\n\n")

############################################################
# PCA Plot
############################################################

cat("Generating PCA Plot...\n")

pcaData <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(
    pcaData,
    aes(PC1, PC2, color = Group)
) +
    geom_point(size = 4) +
    xlab(paste0("PC1 (", percentVar[1], "%)")) +
    ylab(paste0("PC2 (", percentVar[2], "%)")) +
    ggtitle("PCA Plot") +
    theme_classic(base_size = 14)

ggsave("results/plots/PCA_plot.png", p, width = 8, height = 6, dpi = 300)

cat("Done.\n\n")

############################################################
# Sample Distance Heatmap
############################################################

cat("Generating Sample Distance Heatmap...\n")

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- metadata_filtered$Sample
colnames(sampleDistMatrix) <- metadata_filtered$Sample

pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    color = colorRampPalette(
        rev(brewer.pal(9, "Blues"))
    )(255),
    filename = "results/plots/Sample_Distance_Heatmap.png"
)

cat("Done.\n\n")

############################################################
# MA Plot
############################################################

cat("Generating MA Plot...\n")

png(
    "results/plots/MA_plot.png",
    width = 1800,
    height = 1500,
    res = 300
)

plotMA(
    res,
    ylim = c(-6, 6)
)

dev.off()

cat("Done.\n\n")

############################################################
# Volcano Plot
############################################################

cat("Generating Volcano Plot...\n")

########### V2 UPGRADE ###########
# --- Configuration ---
TOP_LABELS <- 15          # Number of top genes to label on the plot
SIG_PADJ   <- 0.05        # Adjusted p-value threshold
SIG_LFC    <- 1           # Absolute log2 fold-change threshold

# --- Prepare volcano data ---
volcano <- res_df
volcano$Gene <- volcano$GeneSymbol
volcano <- volcano[!is.na(volcano$padj), ]

# --- Classify significance ---
volcano$Significance <- "Not Significant"
volcano$Significance[
    volcano$padj < SIG_PADJ & volcano$log2FoldChange >  SIG_LFC
] <- "Upregulated"
volcano$Significance[
    volcano$padj < SIG_PADJ & volcano$log2FoldChange < -SIG_LFC
] <- "Downregulated"

# --- Select top genes for labeling ---
# Rank by significance (smallest padj first), then by absolute fold change
volcano$rank <- order(volcano$padj, decreasing = FALSE)
top_genes <- head(volcano[order(volcano$padj), ], TOP_LABELS)
top_genes <- top_genes[!is.na(top_genes$GeneSymbol) & top_genes$GeneSymbol != "", ]

# --- Volcano summary ---
n_up   <- sum(volcano$Significance == "Upregulated")
n_down <- sum(volcano$Significance == "Downregulated")
n_ns   <- sum(volcano$Significance == "Not Significant")

cat("\n--- Volcano Plot Summary ---\n")
cat("Total genes plotted :", nrow(volcano), "\n")
cat("Upregulated       :", n_up, "\n")
cat("Downregulated     :", n_down, "\n")
cat("Not Significant   :", n_ns, "\n")
cat("Genes labeled     :", nrow(top_genes), "\n")
cat("Thresholds        : padj <", SIG_PADJ, ", |log2FC| >", SIG_LFC, "\n")
cat("-----------------------------\n\n")

# --- Publication-quality volcano plot ---
p <- ggplot(
    volcano,
    aes(
        x = log2FoldChange,
        y = -log10(padj),
        color = Significance
    )
) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_vline(
        xintercept = c(-SIG_LFC, SIG_LFC),
        linetype = "dashed",
        color = "grey40",
        linewidth = 0.6
    ) +
    geom_hline(
        yintercept = -log10(SIG_PADJ),
        linetype = "dashed",
        color = "grey40",
        linewidth = 0.6
    ) +
    geom_text_repel(
        data = top_genes,
        aes(label = Gene),
        size          = 3,
        fontface      = "bold",
        max.overlaps  = Inf,
        box.padding   = 0.4,
        point.padding = 0.3,
        segment.color = "grey50",
        segment.size  = 0.3,
        arrow         = arrow(length = unit(0.01, "npc")),
        show.legend   = FALSE
    ) +
    scale_color_manual(
        name   = "Significance",
        values = c(
            "Downregulated"     = "#2166AC",
            "Not Significant"   = "grey70",
            "Upregulated"     = "#B2182B"
        )
    ) +
    labs(
        title    = "Volcano Plot of Differential Expression",
        subtitle = paste0(
            "n = ", nrow(volcano),
            " | Up: ", n_up,
            " | Down: ", n_down
        ),
        x = expression(log[2] ~ "Fold Change"),
        y = expression(-log[10] ~ "Adjusted P-value")
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        legend.position  = "right",
        legend.title     = element_text(face = "bold"),
        panel.grid.minor = element_blank()
    )

# --- Export high-resolution PNG ---
ggsave("results/plots/Volcano_plot.png", p, width = 8, height = 6, dpi = 300)

# --- Export PDF ---
ggsave("results/plots/Volcano_plot.pdf", p, width = 8, height = 6, device = cairo_pdf)

# --- Export labeled genes CSV ---
write.csv(
    top_genes[, c("GeneSymbol", "GeneID", "log2FoldChange", "pvalue", "padj", "Significance")],
    "results/tables/Volcano_Labeled_Genes.csv",
    row.names = FALSE
)

######## END V2 UPGRADE ########

cat("Done.\n\n")

############################################################
# Heatmap of Top 50 Differentially Expressed Genes
############################################################

cat("Generating Heatmap...\n")

top50 <- head(order(res$padj), 50)

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

cat("Done.\n\n")

############################################################
# Export Gene List for DAVID
############################################################

cat("Exporting DAVID Gene List...\n")

write.table(
    na.omit(sig_genes$GeneSymbol),
    "results/tables/DAVID_Gene_List.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

cat("Done.\n\n")

############################################################
# Summary
############################################################

cat("========================================\n")
cat(" Pipeline Completed Successfully\n")
cat("========================================\n\n")

summary(res)

cat("\nResults saved in:\n")
cat("results/tables/\n")
cat("results/plots/\n")

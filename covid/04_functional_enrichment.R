# 04_functional_enrichment.R
# Functional Enrichment Analysis Pipeline
# GSE201530 RNA-seq
#
# INPUT:
#   results/DESeq2_Results.csv
#
# OUTPUT:
#   GO (BP/MF/CC)
#   KEGG
#   Dotplots
#   Barplots
#   Enrichment maps
#   Gene-concept networks
#   Treeplot
#   Pathview pathway diagram
#
# NOTE:
# Expected columns in DESeq2_Results.csv:
# Gene, log2FoldChange, padj
#
# See clusterProfiler documentation for package installation if needed.

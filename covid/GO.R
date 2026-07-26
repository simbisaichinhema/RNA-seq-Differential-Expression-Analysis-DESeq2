###############################################################
# 04_functional_enrichment.R
# Functional Enrichment Analysis Pipeline
# Project: GSE201530 RNA-seq Analysis
# Author: Simbisai Chinhema
###############################################################

#############################
# 1. Install Required Packages
#############################

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

cran_packages <- c(
    "ggplot2",
    "dplyr",
    "readr"
)

bioc_packages <- c(
    "clusterProfiler",
    "org.Hs.eg.db",
    "AnnotationDbi",
    "enrichplot",
    "DOSE",
    "pathview"
)

for(pkg in cran_packages){

    if(!requireNamespace(pkg, quietly = TRUE)){
        install.packages(pkg)
    }

    library(pkg, character.only = TRUE)
}

for(pkg in bioc_packages){

    if(!requireNamespace(pkg, quietly = TRUE)){
        BiocManager::install(pkg,
                             ask = FALSE,
                             update = FALSE)
    }

    library(pkg, character.only = TRUE)
}

#############################
# 2. Create Output Directories
#############################

dir.create("results", showWarnings = FALSE)

dir.create("results/go",
           recursive = TRUE,
           showWarnings = FALSE)

dir.create("results/kegg",
           recursive = TRUE,
           showWarnings = FALSE)

dir.create("results/plots",
           recursive = TRUE,
           showWarnings = FALSE)

cat("Folders created successfully.\n")

#############################
# 3. Read DESeq2 Results
#############################

cat("Reading DESeq2 results...\n")

res <- read.csv(
    "DESeq2_Results.csv",
    stringsAsFactors = FALSE
)

required_columns <- c(
    "GeneSymbol",
    "GeneID",
    "log2FoldChange",
    "padj"
)

missing_columns <- setdiff(
    required_columns,
    colnames(res)
)

if(length(missing_columns) > 0){

    stop(
        paste(
            "Missing columns:",
            paste(missing_columns,
                  collapse = ", ")
        )
    )
}

res <- res %>%
    filter(
        !is.na(GeneID),
        !is.na(log2FoldChange),
        !is.na(padj)
    )

cat("Total genes loaded:",
    nrow(res),
    "\n")
#############################
# 4. Filter Significant Genes
#############################

cat("Filtering significant genes...\n")

upregulated <- res %>%
    filter(
        padj < 0.05,
        log2FoldChange > 1
    )

downregulated <- res %>%
    filter(
        padj < 0.05,
        log2FoldChange < -1
    )

write.csv(
    upregulated,
    "results/go/Upregulated_Genes.csv",
    row.names = FALSE
)

write.csv(
    downregulated,
    "results/go/Downregulated_Genes.csv",
    row.names = FALSE
)

cat("Upregulated genes :", nrow(upregulated), "\n")
cat("Downregulated genes :", nrow(downregulated), "\n")

#############################
# 5. Prepare Gene List
#############################

gene_list <- unique(
    as.character(upregulated$GeneID)
)

#############################
# 6. GO Biological Process
#############################

cat("Running GO Biological Process...\n")

go_bp <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
)

#############################
# 7. GO Molecular Function
#############################

cat("Running GO Molecular Function...\n")

go_mf <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
)

#############################
# 8. GO Cellular Component
#############################

cat("Running GO Cellular Component...\n")

go_cc <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
)

#############################
# 9. KEGG Pathway Enrichment
#############################

cat("Running KEGG enrichment...\n")

kegg <- enrichKEGG(
    gene          = gene_list,
    organism      = "hsa",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH"
)

#############################
# 10. Save Enrichment Tables
#############################

write.csv(
    as.data.frame(go_bp),
    "results/go/GO_BP.csv",
    row.names = FALSE
)

write.csv(
    as.data.frame(go_mf),
    "results/go/GO_MF.csv",
    row.names = FALSE
)

write.csv(
    as.data.frame(go_cc),
    "results/go/GO_CC.csv",
    row.names = FALSE
)

write.csv(
    as.data.frame(kegg),
    "results/kegg/KEGG.csv",
    row.names = FALSE
)

cat("Enrichment tables saved.\n")

#############################
# 11. Plot Saving Function
#############################

save_plot <- function(filename, plot_object){

    png(
        filename = filename,
        width = 2400,
        height = 1800,
        res = 300
    )

    print(plot_object)

    dev.off()

    pdf(
        file = sub("\\.png$", ".pdf", filename),
        width = 12,
        height = 9
    )

    print(plot_object)

    dev.off()
}

#############################
# 12. GO Dotplots
#############################

cat("Creating GO dotplots...\n")

save_plot(
    "results/plots/GO_BP_Dotplot.png",
    dotplot(
        go_bp,
        showCategory = 20,
        title = "GO Biological Process"
    )
)

save_plot(
    "results/plots/GO_MF_Dotplot.png",
    dotplot(
        go_mf,
        showCategory = 20,
        title = "GO Molecular Function"
    )
)

save_plot(
    "results/plots/GO_CC_Dotplot.png",
    dotplot(
        go_cc,
        showCategory = 20,
        title = "GO Cellular Component"
    )
)

#############################
# 13. GO Barplots
#############################

save_plot(
    "results/plots/GO_BP_Barplot.png",
    barplot(
        go_bp,
        showCategory = 20,
        title = "GO Biological Process"
    )
)

save_plot(
    "results/plots/GO_MF_Barplot.png",
    barplot(
        go_mf,
        showCategory = 20,
        title = "GO Molecular Function"
    )
)

save_plot(
    "results/plots/GO_CC_Barplot.png",
    barplot(
        go_cc,
        showCategory = 20,
        title = "GO Cellular Component"
    )
)

#############################
# 14. KEGG Plots
#############################

cat("Creating KEGG plots...\n")

save_plot(
    "results/plots/KEGG_Dotplot.png",
    dotplot(
        kegg,
        showCategory = 20,
        title = "KEGG Pathway Enrichment"
    )
)

save_plot(
    "results/plots/KEGG_Barplot.png",
    barplot(
        kegg,
        showCategory = 20,
        title = "KEGG Pathway Enrichment"
    )
)

#############################
# 15. GO Enrichment Map
#############################

cat("Creating GO enrichment network plots...\n")

bp_sim <- tryCatch(
    pairwise_termsim(go_bp),
    error = function(e) NULL
)

if(!is.null(bp_sim)){

    save_plot(
        "results/plots/GO_BP_EnrichmentMap.png",
        emapplot(
            bp_sim,
            showCategory = 20
        )
    )

  save_plot(
    "results/plots/GO_BP_Cnetplot.png",
    cnetplot(
      bp_sim,
      showCategory = 10
    )
  )

    save_plot(
        "results/plots/GO_BP_Treeplot.png",
        treeplot(bp_sim)
    )

}else{

    cat("GO enrichment map skipped.\n")

}

#############################
# 16. KEGG Network Plots
#############################

cat("Creating KEGG network plots...\n")

kegg_sim <- tryCatch(
    pairwise_termsim(kegg),
    error = function(e) NULL
)

if(!is.null(kegg_sim)){

    save_plot(
        "results/plots/KEGG_EnrichmentMap.png",
        emapplot(
            kegg_sim,
            showCategory = 20
        )
    )

  save_plot(
    "results/plots/KEGG_Cnetplot.png",
    cnetplot(
      kegg_sim,
      showCategory = 10
    )
    )

}else{

    cat("KEGG enrichment map skipped.\n")

}

#############################
# 17. Pathview
#############################

cat("Generating KEGG pathway diagram...\n")

fold_change <- upregulated$log2FoldChange
names(fold_change) <- as.character(upregulated$GeneID)

if(nrow(as.data.frame(kegg)) > 0){

    top_pathway <- as.data.frame(kegg)$ID[1]

    tryCatch(

        pathview(
            gene.data = fold_change,
            pathway.id = top_pathway,
            species = "hsa",
            out.suffix = "RNAseq"
        ),

        error = function(e){

            cat("Pathview failed:\n")
            print(e$message)

        }

    )

}else{

    cat("No enriched KEGG pathway found.\n")

}

#############################
# 18. Save Session Information
#############################

capture.output(

    sessionInfo(),

    file = "results/sessionInfo.txt"

)

#############################
# 19. Create Analysis Summary
#############################

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

write.csv(

    summary_report,

    "results/Analysis_Summary.csv",

    row.names = FALSE

)

#############################
# 20. Finished
#############################

cat("\n========================================\n")
cat(" RNA-seq Functional Enrichment Complete\n")
cat("========================================\n")

cat("Total genes analysed      :", nrow(res), "\n")
cat("Upregulated genes         :", nrow(upregulated), "\n")
cat("Downregulated genes       :", nrow(downregulated), "\n")
cat("GO BP terms              :", nrow(as.data.frame(go_bp)), "\n")
cat("GO MF terms              :", nrow(as.data.frame(go_mf)), "\n")
cat("GO CC terms              :", nrow(as.data.frame(go_cc)), "\n")
cat("KEGG pathways            :", nrow(as.data.frame(kegg)), "\n")

cat("\nAll results have been saved inside the 'results' folder.\n")
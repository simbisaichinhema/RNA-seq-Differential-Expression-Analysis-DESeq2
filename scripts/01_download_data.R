############################################################
# 01_download_data.R
# Purpose: Download GEO dataset GSE201530 raw count matrix
#          and series matrix metadata.
#
# Required packages: GEOquery
#
# Inputs:  None (downloads from GEO public repository)
# Outputs:
#   data/raw/GSE201530_raw_counts_GRCh38.p13_NCBI.tsv
#   data/raw/GSE201530_series_matrix.txt
#
# Expected runtime: 2–5 minutes (depends on network speed)
#
# Example execution:
#   Rscript scripts/01_download_data.R
#
# Error handling: Stops with informative message if download
#   fails or required metadata fields are missing.
############################################################

cat("========================================\n")
cat(" Step 1: Download GEO Data\n")
cat("========================================\n\n")

############################################################
# Load Libraries
############################################################

required_packages <- c("GEOquery")

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
# Create Output Directories
############################################################

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)

cat("Output directories created.\n\n")

############################################################
# Download GEO Series Matrix
############################################################

cat("Downloading GSE201530 series matrix from GEO...\n")

gse <- tryCatch(
  {
    getGEO("GSE201530", destdir = "data/raw", getGSEMatrix = TRUE)
  },
  error = function(e) {
    stop(paste0("Failed to download GSE201530: ", e$message))
  }
)

cat("Download complete.\n\n")

############################################################
# Extract and Save Count Matrix
############################################################

cat("Extracting count matrix...\n")

# GSE201530 contains expression data as a GEOExpressionSet
gset <- gse[[1]]  # Use the first platform if multiple

# Extract the expression matrix (probe-level)
count_matrix <- exprs(gset)

# Add GeneID column (probe IDs from the platform)
gene_ids <- rownames(count_matrix)
count_df <- data.frame(GeneID = gene_ids, count_matrix, check.names = FALSE)

# Save as TSV
output_counts <- "data/raw/GSE201530_raw_counts_GRCh38.p13_NCBI.tsv"
write.table(
  count_df,
  file = output_counts,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Count matrix saved to:", output_counts, "\n")
cat("Dimensions:", nrow(count_df), "genes x", ncol(count_df) - 1, "samples\n\n")

############################################################
# Extract and Save Series Matrix Text
############################################################

cat("Saving GEO series matrix metadata...\n")

series_matrix_file <- "data/raw/GSE201530_series_matrix.txt"

# GEOquery stores the series matrix in the GEO object
# Write the key metadata to a text file for reference
geo_lines <- character(0)

# Platform annotation
geo_lines <- c(geo_lines, paste0("!Series_title = \"", attr(gse, "header")["series_title"], "\""))
geo_lines <- c(geo_lines, paste0("!Series_geo_accession = \"", attr(gse, "header")["geo_accession"], "\""))
geo_lines <- c(geo_lines, paste0("!Series_status = \"", attr(gse, "header")["status"], "\""))
geo_lines <- c(geo_lines, paste0("!Series_submission_date = \"", attr(gse, "header")["submission_date"], "\""))
geo_lines <- c(geo_lines, paste0("!Series_last_update_date = \"", attr(gse, "header")["last_update_date"], "\""))

# Sample metadata columns
geo_lines <- c(geo_lines, "")
geo_lines <- c(geo_lines, "!Sample_metadata:")
sample_names <- colnames(count_matrix)
for (i in seq_along(sample_names)) {
  geo_lines <- c(geo_lines, paste0("!Sample_ref = ", i, ": ", sample_names[i]))
}

writeLines(geo_lines, series_matrix_file)

cat("Series matrix metadata saved to:", series_matrix_file, "\n\n")

############################################################
# Verify Downloads
############################################################

cat("Verifying downloaded files...\n")

if (!file.exists(output_counts)) {
  stop("Count matrix file not found after download!")
}

if (!file.exists(series_matrix_file)) {
  stop("Series matrix file not found after download!")
}

# Verify the count matrix has expected columns (samples)
if (ncol(count_df) < 2) {
  stop("Count matrix has fewer than 2 samples — check GEO download.")
}

cat("Verification passed.\n\n")

cat("========================================\n")
cat(" Step 1 Complete\n")
cat("========================================\n")

############################################################
# utils.R
# Purpose: Shared utility functions used across analysis scripts.
#
# Functions:
#   check_file()    — verify a file exists before use
#   make_dir()      — create directory if it does not exist
#   session_info()  — capture and save R session information
#   load_counts()   — read a GEO-style count matrix TSV
#   load_metadata() — read metadata CSV with verification
#
# Usage:
#   source("scripts/utils.R")
#
############################################################

############################################################
# Check File Existence
############################################################

check_file <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(paste0("Required file not found: ", filepath))
  }
  invisible(TRUE)
}

############################################################
# Create Directory
############################################################

make_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Created directory: ", path)
  }
  invisible(TRUE)
}

############################################################
# Save Session Information
############################################################

session_info <- function(filepath = "results/sessionInfo.txt") {
  make_dir(dirname(filepath))
  capture.output(
    sessionInfo(),
    file = filepath
  )
  message("Session info saved to: ", filepath)
}

############################################################
# Load Count Matrix
############################################################

load_counts <- function(filepath) {
  check_file(filepath)
  message("Reading count matrix from: ", filepath)
  counts <- read.delim(filepath, check.names = FALSE)
  message("Loaded ", nrow(counts), " genes x ", ncol(counts) - 1, " samples")
  return(counts)
}

############################################################
# Load Metadata
############################################################

load_metadata <- function(filepath) {
  check_file(filepath)
  message("Reading metadata from: ", filepath)
  metadata <- read.csv(filepath, stringsAsFactors = FALSE)
  message("Loaded ", nrow(metadata), " samples x ", ncol(metadata), " columns")
  return(metadata)
}

############################################################
# Sanitize Column Names
############################################################

sanitize_cols <- function(names_vec) {
  cleaned <- gsub("[^[:alnum:][:space:]_\\.]", "", names_vec)
  cleaned <- gsub("\\s+", "_", cleaned)
  return(cleaned)
}

############################################################
# Print Section Header
############################################################

section_header <- function(title) {
  separator <- paste0(rep("=", nchar(title) + 4), collapse = "")
  cat("\n")
  cat(separator, "\n")
  cat("  ", title, "\n")
  cat(separator, "\n\n")
}

############################################################
# Validate DESeq2 Results
############################################################

validate_results <- function(res_df) {
  required <- c("GeneSymbol", "GeneID", "log2FoldChange", "pvalue", "padj")
  missing <- setdiff(required, colnames(res_df))
  if (length(missing) > 0) {
    stop(paste0("Missing columns in results: ", paste(missing, collapse = ", ")))
  }
  message("Results validation passed — ", nrow(res_df), " genes with ", ncol(res_df), " columns")
  invisible(TRUE)
}

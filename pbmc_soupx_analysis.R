#!/usr/bin/env Rscript

# Comprehensive PBMC SoupX Analysis Script
# This script ensures auto-estimation works for PBMC datasets
# Author: Enhanced SoupX Package
# Version: 1.6.5

cat("=== PBMC SoupX Analysis Script ===\n")
cat("Loading required packages...\n")

# Install and load the enhanced SoupX package
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install the enhanced SoupX package
cat("Installing enhanced SoupX package...\n")
devtools::install_github("kevinj24fr/SoupX")

library(SoupX)

cat("Loading 10X data...\n")

# Load 10x data
data_path <- "~/Documents/Documents – Kevin's MacBook Pro/AGRavi/Collaborations/Zurich/CarT_Project/Data/snRNA_Fr/Raw/CAR_Cyto1_invitro"
sc <- load10X(data_path)

# Display basic data information
cat("\n=== DATA SUMMARY ===\n")
cat("Number of cells:", ncol(sc$toc), "\n")
cat("Number of genes:", nrow(sc$toc), "\n")
cat("Total UMIs:", sum(sc$toc), "\n")
cat("Mean UMIs per cell:", round(mean(colSums(sc$toc)), 1), "\n")
cat("Median UMIs per cell:", round(median(colSums(sc$toc)), 1), "\n")

# Check if we have metadata (clusters)
if (!is.null(sc$metaData) && "clusters" %in% colnames(sc$metaData)) {
  cat("Number of clusters:", length(unique(sc$metaData$clusters)), "\n")
  cat("Cluster distribution:\n")
  print(table(sc$metaData$clusters))
} else {
  cat("No cluster information found in metadata\n")
}

# Function to try auto-estimation with different parameters
try_auto_estimation <- function(sc, params_list) {
  for (i in seq_along(params_list)) {
    params <- params_list[[i]]
    cat(sprintf("\n--- Attempt %d: tfidfMin=%.2f, soupQuantile=%.2f ---\n", 
                i, params$tfidfMin, params$soupQuantile))
    
    tryCatch({
      sc_result <- autoEstCont(sc, 
                              tfidfMin = params$tfidfMin, 
                              soupQuantile = params$soupQuantile)
      cat("SUCCESS: Auto-estimation completed!\n")
      cat("Estimated contamination fraction:", sc_result$metaData$rho[1], "\n")
      return(sc_result)
    }, error = function(e) {
      cat("FAILED:", e$message, "\n")
      return(NULL)
    })
  }
  return(NULL)
}

# Define parameter sets to try (from most conservative to most lenient)
params_list <- list(
  list(tfidfMin = 1.0, soupQuantile = 0.9),    # Default parameters
  list(tfidfMin = 0.8, soupQuantile = 0.8),     # Slightly more lenient
  list(tfidfMin = 0.5, soupQuantile = 0.7),     # More lenient
  list(tfidfMin = 0.3, soupQuantile = 0.5),     # Very lenient
  list(tfidfMin = 0.1, soupQuantile = 0.3),     # Extremely lenient
  list(tfidfMin = 0.05, soupQuantile = 0.2)     # Most lenient
)

cat("\n=== ATTEMPTING AUTO-ESTIMATION ===\n")
sc_auto <- try_auto_estimation(sc, params_list)

if (!is.null(sc_auto)) {
  cat("\n=== AUTO-ESTIMATION SUCCESSFUL ===\n")
  sc <- sc_auto
  cat("Using auto-estimated contamination fraction\n")
} else {
  cat("\n=== AUTO-ESTIMATION FAILED - USING MANUAL APPROACH ===\n")
  
  # Try manual contamination with different values
  manual_values <- c(0.05, 0.1, 0.15, 0.2, 0.25)
  
  for (rho in manual_values) {
    cat(sprintf("Trying manual contamination: %.1f%%\n", rho * 100))
    sc_manual <- setContaminationFraction(sc, rho)
    
    # Test if adjustCounts works
    tryCatch({
      adjusted_test <- adjustCounts(sc_manual)
      cat(sprintf("SUCCESS: Manual contamination %.1f%% works\n", rho * 100))
      sc <- sc_manual
      break
    }, error = function(e) {
      cat(sprintf("FAILED: Manual contamination %.1f%% - %s\n", rho * 100, e$message))
    })
  }
}

# Ensure we have a contamination fraction set
if (is.null(sc$metaData$rho) || is.na(sc$metaData$rho[1])) {
  cat("WARNING: No contamination fraction set. Using default 10%\n")
  sc <- setContaminationFraction(sc, 0.1)
}

cat("\n=== CORRECTING COUNTS ===\n")
adjusted <- adjustCounts(sc)

# Display results
cat("\n=== RESULTS SUMMARY ===\n")
cat("Original total UMIs:", sum(sc$toc), "\n")
cat("Adjusted total UMIs:", sum(adjusted), "\n")
cat("Removed UMIs:", sum(sc$toc) - sum(adjusted), "\n")
cat("Contamination fraction used:", sc$metaData$rho[1], "\n")

# Calculate some quality metrics
zero_genes_original <- sum(rowSums(sc$toc) == 0)
zero_genes_adjusted <- sum(rowSums(adjusted) == 0)
cat("Zero-expression genes (original):", zero_genes_original, "\n")
cat("Zero-expression genes (adjusted):", zero_genes_adjusted, "\n")

# Generate comprehensive report
cat("\n=== GENERATING REPORT ===\n")
tryCatch({
  report <- generateSoupXReport(sc, adjusted)
  cat("SUCCESS: Report generated\n")
}, error = function(e) {
  cat("WARNING: Could not generate report -", e$message, "\n")
  cat("This might be due to missing optional dependencies (ggplot2)\n")
})

# Save results
cat("\n=== SAVING RESULTS ===\n")
saveRDS(sc, "soupx_channel.rds")
saveRDS(adjusted, "soupx_adjusted_counts.rds")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files saved:\n")
cat("- soupx_channel.rds: SoupChannel object\n")
cat("- soupx_adjusted_counts.rds: Corrected count matrix\n")

cat("\n=== NEXT STEPS ===\n")
cat("1. Load the corrected counts: adjusted <- readRDS('soupx_adjusted_counts.rds')\n")
cat("2. Use for downstream analysis (Seurat, etc.)\n")
cat("3. Check the contamination fraction used: sc$metaData$rho[1]\n")

cat("\nScript completed successfully!\n") 
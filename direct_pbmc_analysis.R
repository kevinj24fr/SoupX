#!/usr/bin/env Rscript

# Direct PBMC SoupX Analysis Script
# Uses the local enhanced SoupX package without external dependencies

cat("=== Direct PBMC SoupX Analysis ===\n")

# Load the local SoupX package directly
cat("Loading local SoupX package...\n")
library(SoupX)

cat("Loading 10X data...\n")

# Load 10x data
sc <- load10X("~/Documents/Documents – Kevin's MacBook Pro/AGRavi/Collaborations/Zurich/CarT_Project/Data/snRNA_Fr/Raw/CAR_Cyto1_invitro")

cat("Data loaded successfully!\n")
cat("Cells:", ncol(sc$toc), "Genes:", nrow(sc$toc), "\n")

# Try auto-estimation with progressively more lenient parameters
cat("\nAttempting auto-estimation...\n")

auto_success <- FALSE
param_sets <- list(
  list(tfidfMin = 1.0, soupQuantile = 0.9),    # Default
  list(tfidfMin = 0.5, soupQuantile = 0.7),     # More lenient
  list(tfidfMin = 0.3, soupQuantile = 0.5),     # Very lenient
  list(tfidfMin = 0.1, soupQuantile = 0.3)      # Extremely lenient
)

for (i in seq_along(param_sets)) {
  params <- param_sets[[i]]
  cat(sprintf("Trying parameters %d: tfidfMin=%.1f, soupQuantile=%.1f\n", 
              i, params$tfidfMin, params$soupQuantile))
  
  tryCatch({
    sc <- autoEstCont(sc, tfidfMin = params$tfidfMin, soupQuantile = params$soupQuantile)
    cat("SUCCESS: Auto-estimation completed!\n")
    cat("Estimated contamination:", sc$metaData$rho[1], "\n")
    auto_success <- TRUE
    break
  }, error = function(e) {
    cat("Failed:", e$message, "\n")
  })
}

# If auto-estimation fails, use manual approach
if (!auto_success) {
  cat("\nAuto-estimation failed. Using manual contamination setting...\n")
  sc <- setContaminationFraction(sc, 0.1)  # 10% contamination
  cat("Manual contamination set to 10%\n")
}

# Correct counts
cat("\nCorrecting counts...\n")
adjusted <- adjustCounts(sc)

cat("Count correction completed!\n")
cat("Original UMIs:", sum(sc$toc), "\n")
cat("Adjusted UMIs:", sum(adjusted), "\n")
cat("Removed UMIs:", sum(sc$toc) - sum(adjusted), "\n")

# Generate comprehensive report
cat("\nGenerating report...\n")
tryCatch({
  report <- generateSoupXReport(sc, adjusted)
  cat("Report generated successfully!\n")
}, error = function(e) {
  cat("Warning: Could not generate report -", e$message, "\n")
})

cat("\nAnalysis completed successfully!\n") 
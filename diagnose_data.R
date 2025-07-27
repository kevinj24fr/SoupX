#!/usr/bin/env Rscript

# Diagnostic script for SoupX data analysis
# This will help understand why autoEstCont is failing

library(SoupX)

# Load the data
cat("Loading 10X data...\n")
sc <- load10X("~/Documents/Documents – Kevin's MacBook Pro/AGRavi/Collaborations/Zurich/CarT_Project/Data/snRNA_Fr/Raw/CAR_Cyto1_invitro")

# Basic data info
cat("\n=== DATA SUMMARY ===\n")
cat("Number of cells:", ncol(sc$toc), "\n")
cat("Number of genes:", nrow(sc$toc), "\n")
cat("Total UMIs:", sum(sc$toc), "\n")
cat("Mean UMIs per cell:", mean(colSums(sc$toc)), "\n")
cat("Median UMIs per cell:", median(colSums(sc$toc)), "\n")

# Check for metadata
if(!is.null(sc$metaData)) {
  cat("Metadata columns:", colnames(sc$metaData), "\n")
  if("clusters" %in% colnames(sc$metaData)) {
    cat("Number of clusters:", length(unique(sc$metaData$clusters)), "\n")
    cat("Cluster sizes:\n")
    print(table(sc$metaData$clusters))
  }
}

# Check soup profile
if(!is.null(sc$soupProfile)) {
  cat("\n=== SOUP PROFILE ===\n")
  cat("Top 10 soup genes:\n")
  soup_genes <- head(sort(sc$soupProfile, decreasing=TRUE), 10)
  print(soup_genes)
}

# Try different parameters for autoEstCont
cat("\n=== TESTING AUTO-ESTIMATION PARAMETERS ===\n")

# Test 1: Very lenient parameters
cat("Test 1: Very lenient parameters (tfidfMin=0.1, soupQuantile=0.5)\n")
tryCatch({
  sc_test1 <- autoEstCont(sc, tfidfMin=0.1, soupQuantile=0.5)
  cat("SUCCESS with lenient parameters!\n")
}, error=function(e) {
  cat("FAILED with lenient parameters:", e$message, "\n")
})

# Test 2: Even more lenient
cat("\nTest 2: Extremely lenient parameters (tfidfMin=0.05, soupQuantile=0.3)\n")
tryCatch({
  sc_test2 <- autoEstCont(sc, tfidfMin=0.05, soupQuantile=0.3)
  cat("SUCCESS with extremely lenient parameters!\n")
}, error=function(e) {
  cat("FAILED with extremely lenient parameters:", e$message, "\n")
})

# Test 3: Manual contamination setting
cat("\nTest 3: Manual contamination setting\n")
sc_manual <- setContaminationFraction(sc, 0.1)
cat("SUCCESS: Manual contamination set to 10%\n")

# Test adjustCounts with manual setting
cat("\nTest 4: Adjusting counts with manual contamination\n")
adjusted <- adjustCounts(sc_manual)
cat("SUCCESS: Counts adjusted!\n")
cat("Original total UMIs:", sum(sc$toc), "\n")
cat("Adjusted total UMIs:", sum(adjusted), "\n")
cat("Removed UMIs:", sum(sc$toc) - sum(adjusted), "\n")

cat("\n=== RECOMMENDATION ===\n")
cat("Since auto-estimation is failing, use manual contamination setting:\n")
cat("sc <- setContaminationFraction(sc, 0.1)  # Try 10% contamination\n")
cat("adjusted <- adjustCounts(sc)\n")
cat("You can try different values: 0.05, 0.1, 0.15, 0.2\n") 
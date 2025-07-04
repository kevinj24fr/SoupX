#!/usr/bin/env Rscript

# Test script for single-cluster edge cases in SoupX
# This tests all functions that could fail with single clusters

library(SoupX)

# Create a minimal test dataset with single cluster
set.seed(42)
n_genes <- 100
n_cells <- 50

# Create toy data with single cluster
toc <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.1)
rownames(toc) <- paste0("Gene", 1:n_genes)
colnames(toc) <- paste0("Cell", 1:n_cells)

# Create soup profile
soup_profile <- data.frame(
  est = runif(n_genes, 0, 0.1),
  counts = rpois(n_genes, 10)
)
rownames(soup_profile) <- rownames(toc)

# Create metadata with single cluster
meta_data <- data.frame(
  nUMIs = rpois(n_cells, 1000),
  rho = runif(n_cells, 0.05, 0.15),
  clusters = rep("Cluster1", n_cells)
)
rownames(meta_data) <- colnames(toc)

# Create SoupChannel object
sc <- list(
  toc = toc,
  tod = toc,  # For simplicity, use same matrix
  soupProfile = soup_profile,
  metaData = meta_data
)
class(sc) <- "SoupChannel"

cat("Testing single cluster edge cases...\n")

# Test 1: autoEstCont with single cluster
cat("1. Testing autoEstCont with single cluster...\n")
tryCatch({
  result <- autoEstCont(sc, verbose = FALSE)
  cat("   ✓ autoEstCont passed\n")
}, error = function(e) {
  cat("   ✗ autoEstCont failed:", e$message, "\n")
})

# Test 2: adjustCounts with single cluster
cat("2. Testing adjustCounts with single cluster...\n")
tryCatch({
  result <- adjustCounts(sc, verbose = 0)
  cat("   ✓ adjustCounts passed\n")
}, error = function(e) {
  cat("   ✗ adjustCounts failed:", e$message, "\n")
})

# Test 3: estimateNonExpressingCells with single cluster
cat("3. Testing estimateNonExpressingCells with single cluster...\n")
tryCatch({
  gene_list <- list(Test = c("Gene1", "Gene2", "Gene3"))
  result <- estimateNonExpressingCells(sc, gene_list)
  cat("   ✓ estimateNonExpressingCells passed\n")
}, error = function(e) {
  cat("   ✗ estimateNonExpressingCells failed:", e$message, "\n")
})

# Test 4: expandClusters with single cluster
cat("4. Testing expandClusters with single cluster...\n")
tryCatch({
  # Create dummy cluster-level counts
  clust_counts <- Matrix::rsparsematrix(n_genes, 1, density = 0.1)
  rownames(clust_counts) <- rownames(toc)
  colnames(clust_counts) <- "Cluster1"
  
  clusters <- setNames(rep("Cluster1", n_cells), colnames(toc))
  cell_weights <- meta_data$nUMIs * meta_data$rho
  
  result <- SoupX::expandClusters(clust_counts, toc, clusters, cell_weights, verbose = 0)
  cat("   ✓ expandClusters passed\n")
}, error = function(e) {
  cat("   ✗ expandClusters failed:", e$message, "\n")
})

# Test 5: quickMarkers with single cluster
cat("5. Testing quickMarkers with single cluster...\n")
tryCatch({
  result <- quickMarkers(toc, rep("Cluster1", n_cells), N = 10)
  cat("   ✓ quickMarkers passed\n")
}, error = function(e) {
  cat("   ✗ quickMarkers failed:", e$message, "\n")
})

# Test 6: Full pipeline with single cluster
cat("6. Testing full pipeline with single cluster...\n")
tryCatch({
  # Set contamination fraction manually
  sc <- setContaminationFraction(sc, 0.1)
  
  # Run adjustment
  result <- adjustCounts(sc, verbose = 0)
  cat("   ✓ Full pipeline passed\n")
}, error = function(e) {
  cat("   ✗ Full pipeline failed:", e$message, "\n")
})

# Test 7: Multiple single clusters (edge case)
cat("7. Testing with multiple single-cell clusters...\n")
tryCatch({
  # Create individual clusters for each cell
  individual_clusters <- setNames(paste0("Cell", 1:n_cells), colnames(toc))
  sc$metaData$clusters <- individual_clusters
  
  result <- autoEstCont(sc, verbose = FALSE)
  cat("   ✓ Multiple single-cell clusters passed\n")
}, error = function(e) {
  cat("   ✗ Multiple single-cell clusters failed:", e$message, "\n")
})

cat("\nAll single cluster edge case tests completed!\n") 
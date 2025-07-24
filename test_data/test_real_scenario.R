#!/usr/bin/env Rscript

# Test script to verify the autoEstCont fix works with real data scenario
# This simulates the exact error the user encountered

library(SoupX)

cat("==== Testing Real Scenario Fix ====\n")

# Create data that mimics the user's real data structure
set.seed(42)
n_genes <- 1000
n_cells <- 500

# Create realistic count matrix
toc <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.1)
rownames(toc) <- paste0("Gene", 1:n_genes)
colnames(toc) <- paste0("Cell", 1:n_cells)

# Create soup profile
soup_profile <- data.frame(
  est = runif(n_genes, 0, 0.1),
  counts = rpois(n_genes, 10)
)
rownames(soup_profile) <- rownames(toc)

# Create metadata with single cluster (exactly like user's scenario)
meta_data <- data.frame(
  nUMIs = rpois(n_cells, 1000),
  rho = runif(n_cells, 0.05, 0.15),
  clusters = rep("Cluster1", n_cells)  # Single cluster for all cells
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

cat("Testing autoEstCont with single cluster (user's exact scenario)...\n")

# Test the exact scenario that was failing
tryCatch({
  result <- autoEstCont(sc, verbose = FALSE, doPlot = FALSE)
  cat("✓ SUCCESS: autoEstCont worked without logical subscript error!\n")
  cat("  - Estimated contamination fraction:", result$metaData$rho[1], "\n")
  cat("  - Number of cells processed:", nrow(result$metaData), "\n")
  cat("  - Number of genes processed:", nrow(result$toc), "\n")
}, error = function(e) {
  cat("✗ FAILED: Still getting error:", e$message, "\n")
  cat("  Error type:", class(e)[1], "\n")
  cat("  This means the fix didn't work.\n")
})

cat("\n==== Test completed ====\n") 
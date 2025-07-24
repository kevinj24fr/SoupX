#!/usr/bin/env Rscript

# Minimal test to isolate quickMarkers dimnames error

library(SoupX)

cat("==== Minimal quickMarkers Test ====\n")

# Create minimal test data
set.seed(42)
n_genes <- 5
n_cells <- 3

toc <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.5)
rownames(toc) <- paste0("Gene", 1:n_genes)
colnames(toc) <- paste0("Cell", 1:n_cells)

clusters <- rep("Cluster1", n_cells)

cat("Data created:\n")
cat("  - toc dims:", dim(toc), "\n")
cat("  - toc class:", class(toc), "\n")
cat("  - clusters length:", length(clusters), "\n")

# Test quickMarkers step by step
cat("\nStep 1: Converting to TsparseMatrix...\n")
tryCatch({
  toc_sparse = as(toc, "TsparseMatrix")
  cat("  - Conversion successful\n")
}, error = function(e) {
  cat("  - Conversion failed:", e$message, "\n")
})

cat("\nStep 2: Finding expressed genes...\n")
tryCatch({
  w = which(toc_sparse@x > 0.9)
  cat("  - Expressed genes found:", length(w), "\n")
}, error = function(e) {
  cat("  - Finding expressed genes failed:", e$message, "\n")
})

cat("\nStep 3: Getting cluster counts...\n")
tryCatch({
  clCnts = table(clusters)
  cat("  - Cluster counts:", clCnts, "\n")
}, error = function(e) {
  cat("  - Getting cluster counts failed:", e$message, "\n")
})

cat("\nStep 4: Creating observed counts matrix...\n")
tryCatch({
  gene_names <- rownames(toc_sparse)
  cluster_names <- names(clCnts)
  cat("  - Gene names length:", length(gene_names), "\n")
  cat("  - Cluster names length:", length(cluster_names), "\n")
  
  nObs <- matrix(0, nrow=length(gene_names), ncol=length(cluster_names))
  cat("  - nObs created, dims:", dim(nObs), "\n")
  
  rownames(nObs) <- gene_names
  cat("  - Rownames set\n")
  
  colnames(nObs) <- cluster_names
  cat("  - Colnames set\n")
  
  cat("  - nObs dimnames: rows=", length(rownames(nObs)), ", cols=", length(colnames(nObs)), "\n")
}, error = function(e) {
  cat("  - Creating observed counts matrix failed:", e$message, "\n")
})

cat("\nStep 5: Filling observed counts...\n")
tryCatch({
  for(i in seq_along(cluster_names)) {
    cluster_name <- cluster_names[i]
    cluster_cells <- which(clusters == cluster_name)
    cat("    - Cluster", cluster_name, "has", length(cluster_cells), "cells\n")
    
    if(length(cluster_cells) > 0) {
      cluster_data <- toc_sparse[, cluster_cells, drop=FALSE]
      cat("    - Cluster data dims:", dim(cluster_data), "\n")
      
      if(length(cluster_data@x) > 0) {
        expr_mask <- cluster_data@x > 0.9
        cat("    - Expression mask length:", length(expr_mask), "\n")
        cat("    - Expression mask sum:", sum(expr_mask), "\n")
        
        if(any(expr_mask)) {
          cluster_data@x <- cluster_data@x[expr_mask]
          cluster_data@i <- cluster_data@i[expr_mask]
          cluster_data@j <- cluster_data@j[expr_mask]
          cat("    - Filtered cluster data, kept", sum(expr_mask), "elements\n")
        }
      }
      
      if(length(cluster_data@x) > 0) {
        gene_counts <- table(factor(gene_names[cluster_data@i + 1], levels=gene_names))
        cat("    - Gene counts calculated, length:", length(gene_counts), "\n")
        
        gene_counts <- as.numeric(gene_counts)
        names(gene_counts) <- gene_names
        cat("    - Gene counts named, length:", length(gene_counts), "\n")
        
        nObs[, i] <- gene_counts[gene_names]
        cat("    - Gene counts assigned to nObs\n")
      }
    }
  }
  cat("  - Filling observed counts successful\n")
}, error = function(e) {
  cat("  - Filling observed counts failed:", e$message, "\n")
})

cat("\nStep 6: Testing full quickMarkers...\n")
tryCatch({
  # Test with a lower FDR threshold to get more markers
  mrks = quickMarkers(toc, clusters, N = Inf, FDR = 0.1)
  cat("✓ SUCCESS: quickMarkers worked!\n")
  cat("  - Markers found:", nrow(mrks), "\n")
}, error = function(e) {
  cat("✗ FAILED: quickMarkers error:", e$message, "\n")
  cat("  - Error type:", class(e)[1], "\n")
  cat("  - Error call:", e$call, "\n")
})

cat("\n==== Test completed ====\n") 
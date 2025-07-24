#!/usr/bin/env Rscript

# Minimal test to isolate quickMarkers dimnames error

library(SoupX)

cat("==== Testing quickMarkers dimnames error ====\n")

# Create minimal test data
set.seed(42)
n_genes <- 10
n_cells <- 5

toc <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3)
rownames(toc) <- paste0("Gene", 1:n_genes)
colnames(toc) <- paste0("Cell", 1:n_cells)

clusters <- rep("Cluster1", n_cells)

cat("Data created:\n")
cat("  - toc dims:", dim(toc), "\n")
cat("  - toc class:", class(toc), "\n")
cat("  - clusters length:", length(clusters), "\n")
cat("  - clusters:", clusters, "\n")

# Test quickMarkers directly
cat("\nTesting quickMarkers...\n")
tryCatch({
  # Test step by step
  cat("  - Converting to TsparseMatrix...\n")
  toc_sparse = as(toc, "TsparseMatrix")
  
  cat("  - Finding expressed genes...\n")
  w = which(toc_sparse@x > 0.9)
  cat("  - Expressed genes found:", length(w), "\n")
  
  cat("  - Getting cluster counts...\n")
  clCnts = table(clusters)
  cat("  - Cluster counts:", clCnts, "\n")
  
  cat("  - Creating observed counts matrix...\n")
  gene_names <- rownames(toc_sparse)
  cluster_names <- names(clCnts)
  nObs <- Matrix::Matrix(0, nrow=length(gene_names), ncol=length(cluster_names), sparse=TRUE)
  cat("  - nObs dims:", dim(nObs), "\n")
  
  cat("  - Setting rownames...\n")
  rownames(nObs) <- gene_names
  cat("  - Setting colnames...\n")
  colnames(nObs) <- cluster_names
  
  cat("  - Filling observed counts...\n")
  for(i in seq_along(cluster_names)) {
    cluster_name <- cluster_names[i]
    cluster_cells <- which(clusters == cluster_name)
    cat("    - Cluster", cluster_name, "has", length(cluster_cells), "cells\n")
    if(length(cluster_cells) > 0) {
      cluster_data <- toc_sparse[, cluster_cells, drop=FALSE]
      cat("    - Cluster data dims before filtering:", dim(cluster_data), "\n")
      cat("    - Cluster data x values:", cluster_data@x, "\n")
      cat("    - Values > 0.9:", sum(cluster_data@x > 0.9), "\n")
      # Filter by expression threshold properly for sparse matrices
      if(length(cluster_data@x) > 0) {
        # Create a logical mask for the non-zero elements
        expr_mask <- cluster_data@x > 0.9
        if(any(expr_mask)) {
          # Keep only the elements that pass the threshold
          cluster_data@x <- cluster_data@x[expr_mask]
          cluster_data@i <- cluster_data@i[expr_mask]
          cluster_data@j <- cluster_data@j[expr_mask]
          cat("    - Filtered cluster data, kept", sum(expr_mask), "elements\n")
        } else {
          # No elements pass threshold, create empty matrix
          cluster_data <- cluster_data[, integer(0), drop=FALSE]
          cat("    - No elements pass threshold\n")
        }
      }
      cat("    - Cluster data dims after filtering:", dim(cluster_data), "\n")
              if(length(cluster_data@x) > 0) {
          cat("    - About to calculate gene counts...\n")
          gene_counts <- table(factor(gene_names[cluster_data@i + 1], levels=gene_names))
          cat("    - Gene counts calculated, length:", length(gene_counts), "\n")
          gene_counts <- as.numeric(gene_counts)
          names(gene_counts) <- gene_names
          cat("    - About to assign to nObs...\n")
          nObs[, i] <- gene_counts
          cat("    - Gene counts set for cluster", cluster_name, "\n")
        }
    }
  }
  
  cat("  - Calculating frequencies...\n")
  cat("    - nObs dims:", dim(nObs), "\n")
  cat("    - nObs class:", class(nObs), "\n")
  cat("    - clCnts:", clCnts, "\n")
  cat("    - colnames(nObs):", colnames(nObs), "\n")
  
  cat("    - About to calculate rowSums...\n")
  # Use Matrix::rowSums for sparse matrices
  if(inherits(nObs, "Matrix")) {
    nTot = Matrix::rowSums(nObs)
  } else {
    nTot = rowSums(nObs)
  }
  # Ensure nTot is a vector for single-cluster case
  if(length(cluster_names) == 1) {
    nTot = as.numeric(nTot)
  }
  cat("    - nTot length:", length(nTot), "\n")
  cat("    - nTot class:", class(nTot), "\n")
  
  cat("    - About to calculate tf...\n")
  cat("    - nObs dims for tf:", dim(nObs), "\n")
  cat("    - clCnts[colnames(nObs)]:", clCnts[colnames(nObs)], "\n")
  # Convert sparse matrix to regular matrix for transpose operations
  nObs_mat = as.matrix(nObs)
  tf = t(t(nObs_mat)/as.integer(clCnts[colnames(nObs)]))
  cat("    - tf dims:", dim(tf), "\n")
  
  cat("    - About to calculate ntf...\n")
  # Handle single-cluster case where nTot is a vector
  if(length(cluster_names) == 1) {
    ntf = t(t(nTot - nObs_mat)/as.integer(ncol(toc_sparse)-clCnts[colnames(nObs)]))
  } else {
    ntf = t(t(nTot - nObs_mat)/as.integer(ncol(toc_sparse)-clCnts[colnames(nObs)]))
  }
  cat("    - ntf dims:", dim(ntf), "\n")
  
  idf = log(ncol(toc_sparse)/nTot)
  cat("    - idf length:", length(idf), "\n")
  
  score = tf*idf
  cat("    - score dims:", dim(score), "\n")
  
  cat("  - Calculating p-values...\n")
  qvals = matrix(0, nrow=nrow(nObs), ncol=ncol(nObs))
  for(i in seq_len(ncol(nObs))) {
    pvals <- phyper(nObs[,i]-1, nTot, ncol(toc_sparse)-nTot, clCnts[colnames(nObs)[i]], lower.tail=FALSE)
    qvals[,i] <- p.adjust(pvals, method='BH')
  }
  colnames(qvals) = colnames(nObs)
  
  cat("  - Calculating second best...\n")
  sndBest = matrix(0, nrow=nrow(tf), ncol=ncol(tf))
  sndBestName = matrix("", nrow=nrow(tf), ncol=ncol(tf))
  
  for(i in seq_len(ncol(tf))) {
    other_cols <- setdiff(seq_len(ncol(tf)), i)
    cat("    - Column", i, "other_cols:", other_cols, "\n")
    if(length(other_cols) > 0) {
      other_data <- tf[, other_cols, drop=FALSE]
      cat("    - Other data dims:", dim(other_data), "\n")
      if(ncol(other_data) > 0 && nrow(other_data) > 0) {
        max_vals <- apply(other_data, 1, function(x) {
          if(all(is.na(x)) || all(x == 0)) return(0) else max(x, na.rm=TRUE)
        })
        max_idx <- apply(other_data, 1, function(x) {
          if(all(is.na(x)) || all(x == 0)) return(1) else which.max(x)
        })
        sndBest[, i] <- max_vals
        sndBestName[, i] <- colnames(tf)[other_cols[max_idx]]
        cat("    - Second best calculated for column", i, "\n")
      }
    }
  }
  
  cat("  - Setting matrix dimnames...\n")
  colnames(sndBest) = colnames(tf)
  colnames(sndBestName) = colnames(tf)
  rownames(sndBestName) = rownames(tf)
  
  cat("  - Getting top markers...\n")
  w = lapply(seq_len(ncol(nObs)), function(e){
    o = order(score[,e], decreasing=TRUE)
    if(sum(qvals[,e]<0.01)>=10){
      o[seq(10)]
    }else{
      o[qvals[o,e]<0.01]
    }
  })
  
  cat("  - Constructing output...\n")
  if(length(unlist(w)) == 0) {
    cat("  - No markers found, returning empty result\n")
    return(data.frame())
  }
  
  ww = cbind(unlist(w,use.names=FALSE),rep(seq_len(ncol(nObs)),lengths(w)))
  
  valid_rows = ww[,1] <= nrow(nObs) & ww[,1] > 0
  valid_cols = ww[,2] <= ncol(nObs) & ww[,2] > 0
  valid_indices = valid_rows & valid_cols
  
  if(sum(valid_indices) == 0) {
    cat("  - No valid indices found, returning empty result\n")
    return(data.frame())
  }
  
  ww = ww[valid_indices, , drop = FALSE]
  
  cat("  - Creating final data frame...\n")
  out = data.frame(gene = rownames(nObs)[ww[,1]],
                   cluster = colnames(nObs)[ww[,2]],
                   geneFrequency = tf[ww],
                   geneFrequencyOutsideCluster = ntf[ww],
                   geneFrequencySecondBest = sndBest[ww],
                   geneFrequencyGlobal = nTot[ww[,1]]/ncol(toc_sparse),
                   secondBestClusterName = sndBestName[ww],
                   tfidf = score[ww],
                   idf = idf[ww[,1]],
                   qval = qvals[ww],
                   stringsAsFactors=FALSE)
  
  cat("✓ SUCCESS: quickMarkers worked!\n")
  cat("  - Markers found:", nrow(out), "\n")
}, error = function(e) {
  cat("✗ FAILED: quickMarkers error:", e$message, "\n")
  cat("  - Error type:", class(e)[1], "\n")
})

cat("\n==== Test completed ====\n") 
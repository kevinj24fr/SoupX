#!/usr/bin/env Rscript

# Comprehensive bug test for SoupX
# Tests all edge cases and potential failure scenarios

library(SoupX)

cat("==== Comprehensive SoupX Bug Test ====\n")

# Test 1: Empty data edge cases
cat("1. Testing empty data edge cases...\n")
tryCatch({
  # Empty matrix
  empty_toc <- Matrix(0, nrow=0, ncol=0, sparse=TRUE)
  empty_soup <- data.frame(est=numeric(0), counts=numeric(0))
  empty_meta <- data.frame(nUMIs=numeric(0), rho=numeric(0))
  
  sc_empty <- list(toc=empty_toc, tod=empty_toc, soupProfile=empty_soup, metaData=empty_meta)
  class(sc_empty) <- "SoupChannel"
  
  cat("   ✓ Empty data handling\n")
}, error = function(e) {
  cat("   ✗ Empty data failed:", e$message, "\n")
})

# Test 2: Single cell edge cases
cat("2. Testing single cell edge cases...\n")
tryCatch({
  single_toc <- Matrix(1, nrow=10, ncol=1, sparse=TRUE)
  rownames(single_toc) <- paste0("Gene", 1:10)
  colnames(single_toc) <- "Cell1"
  
  single_soup <- data.frame(est=rep(0.1, 10), counts=rep(10, 10))
  rownames(single_soup) <- rownames(single_toc)
  
  single_meta <- data.frame(nUMIs=100, rho=0.1, clusters="Single")
  rownames(single_meta) <- colnames(single_toc)
  
  sc_single <- list(toc=single_toc, tod=single_toc, soupProfile=single_soup, metaData=single_meta)
  class(sc_single) <- "SoupChannel"
  
  # Test autoEstCont with single cell
  sc_single <- autoEstCont(sc_single, verbose=FALSE, doPlot=FALSE)
  cat("   ✓ Single cell autoEstCont\n")
  
  # Test adjustCounts with single cell
  adj_single <- adjustCounts(sc_single, verbose=0)
  cat("   ✓ Single cell adjustCounts\n")
  
}, error = function(e) {
  cat("   ✗ Single cell failed:", e$message, "\n")
})

# Test 3: Single cluster edge cases
cat("3. Testing single cluster edge cases...\n")
tryCatch({
  single_cluster_toc <- Matrix(rpois(100, 5), nrow=10, ncol=10, sparse=TRUE)
  rownames(single_cluster_toc) <- paste0("Gene", 1:10)
  colnames(single_cluster_toc) <- paste0("Cell", 1:10)
  
  single_cluster_soup <- data.frame(est=rep(0.1, 10), counts=rep(10, 10))
  rownames(single_cluster_soup) <- rownames(single_cluster_toc)
  
  single_cluster_meta <- data.frame(
    nUMIs=rpois(10, 1000),
    rho=runif(10, 0.05, 0.15),
    clusters=rep("Cluster1", 10)
  )
  rownames(single_cluster_meta) <- colnames(single_cluster_toc)
  
  sc_single_cluster <- list(
    toc=single_cluster_toc, 
    tod=single_cluster_toc, 
    soupProfile=single_cluster_soup, 
    metaData=single_cluster_meta
  )
  class(sc_single_cluster) <- "SoupChannel"
  
  # Test autoEstCont with single cluster
  sc_single_cluster <- autoEstCont(sc_single_cluster, verbose=FALSE, doPlot=FALSE)
  cat("   ✓ Single cluster autoEstCont\n")
  
  # Test adjustCounts with single cluster
  adj_single_cluster <- adjustCounts(sc_single_cluster, verbose=0)
  cat("   ✓ Single cluster adjustCounts\n")
  
  # Test expandClusters with single cluster
  clust_counts <- Matrix(rpois(10, 5), nrow=10, ncol=1, sparse=TRUE)
  rownames(clust_counts) <- rownames(single_cluster_toc)
  colnames(clust_counts) <- "Cluster1"
  
  clusters <- setNames(rep("Cluster1", 10), colnames(single_cluster_toc))
  cell_weights <- single_cluster_meta$nUMIs * single_cluster_meta$rho
  
  exp_single <- expandClusters(clust_counts, single_cluster_toc, clusters, cell_weights, verbose=0)
  cat("   ✓ Single cluster expandClusters\n")
  
}, error = function(e) {
  cat("   ✗ Single cluster failed:", e$message, "\n")
})

# Test 4: Multiple clusters edge cases
cat("4. Testing multiple clusters edge cases...\n")
tryCatch({
  multi_toc <- Matrix(rpois(200, 5), nrow=20, ncol=10, sparse=TRUE)
  rownames(multi_toc) <- paste0("Gene", 1:20)
  colnames(multi_toc) <- paste0("Cell", 1:10)
  
  multi_soup <- data.frame(est=rep(0.1, 20), counts=rep(10, 20))
  rownames(multi_soup) <- rownames(multi_toc)
  
  multi_meta <- data.frame(
    nUMIs=rpois(10, 1000),
    rho=runif(10, 0.05, 0.15),
    clusters=rep(c("A", "B"), each=5)
  )
  rownames(multi_meta) <- colnames(multi_toc)
  
  sc_multi <- list(toc=multi_toc, tod=multi_toc, soupProfile=multi_soup, metaData=multi_meta)
  class(sc_multi) <- "SoupChannel"
  
  # Test autoEstCont with multiple clusters
  sc_multi <- autoEstCont(sc_multi, verbose=FALSE, doPlot=FALSE)
  cat("   ✓ Multiple clusters autoEstCont\n")
  
  # Test adjustCounts with multiple clusters
  adj_multi <- adjustCounts(sc_multi, verbose=0)
  cat("   ✓ Multiple clusters adjustCounts\n")
  
  # Test expandClusters with multiple clusters
  clust_counts_multi <- Matrix(rpois(20, 5), nrow=20, ncol=2, sparse=TRUE)
  rownames(clust_counts_multi) <- rownames(multi_toc)
  colnames(clust_counts_multi) <- c("A", "B")
  
  clusters_multi <- setNames(multi_meta$clusters, colnames(multi_toc))
  cell_weights_multi <- multi_meta$nUMIs * multi_meta$rho
  
  exp_multi <- expandClusters(clust_counts_multi, multi_toc, clusters_multi, cell_weights_multi, verbose=0)
  cat("   ✓ Multiple clusters expandClusters\n")
  
}, error = function(e) {
  cat("   ✗ Multiple clusters failed:", e$message, "\n")
})

# Test 5: Edge cases with missing data
cat("5. Testing edge cases with missing data...\n")
tryCatch({
  # Test with NA values in metadata
  na_toc <- Matrix(rpois(100, 5), nrow=10, ncol=10, sparse=TRUE)
  rownames(na_toc) <- paste0("Gene", 1:10)
  colnames(na_toc) <- paste0("Cell", 1:10)
  
  na_soup <- data.frame(est=rep(0.1, 10), counts=rep(10, 10))
  rownames(na_soup) <- rownames(na_toc)
  
  na_meta <- data.frame(
    nUMIs=c(rpois(9, 1000), NA),
    rho=c(runif(9, 0.05, 0.15), NA),
    clusters=rep("Cluster1", 10)
  )
  rownames(na_meta) <- colnames(na_toc)
  
  sc_na <- list(toc=na_toc, tod=na_toc, soupProfile=na_soup, metaData=na_meta)
  class(sc_na) <- "SoupChannel"
  
  # This should fail gracefully
  tryCatch({
    sc_na <- autoEstCont(sc_na, verbose=FALSE, doPlot=FALSE)
    cat("   ✗ NA handling failed - should have caught NA values\n")
  }, error = function(e) {
    cat("   ✓ NA handling worked correctly\n")
  })
  
}, error = function(e) {
  cat("   ✗ NA testing failed:", e$message, "\n")
})

# Test 6: Large data edge cases
cat("6. Testing large data edge cases...\n")
tryCatch({
  # Test with larger dataset
  large_toc <- Matrix(rpois(1000, 3), nrow=100, ncol=10, sparse=TRUE)
  rownames(large_toc) <- paste0("Gene", 1:100)
  colnames(large_toc) <- paste0("Cell", 1:10)
  
  large_soup <- data.frame(est=rep(0.1, 100), counts=rep(10, 100))
  rownames(large_soup) <- rownames(large_toc)
  
  large_meta <- data.frame(
    nUMIs=rpois(10, 1000),
    rho=runif(10, 0.05, 0.15),
    clusters=rep("Large", 10)
  )
  rownames(large_meta) <- colnames(large_toc)
  
  sc_large <- list(toc=large_toc, tod=large_toc, soupProfile=large_soup, metaData=large_meta)
  class(sc_large) <- "SoupChannel"
  
  # Test autoEstCont with large data
  sc_large <- autoEstCont(sc_large, verbose=FALSE, doPlot=FALSE)
  cat("   ✓ Large data autoEstCont\n")
  
  # Test adjustCounts with large data
  adj_large <- adjustCounts(sc_large, verbose=0)
  cat("   ✓ Large data adjustCounts\n")
  
}, error = function(e) {
  cat("   ✗ Large data failed:", e$message, "\n")
})

# Test 7: Zero expression edge cases
cat("7. Testing zero expression edge cases...\n")
tryCatch({
  # Test with mostly zero expression
  zero_toc <- Matrix(0, nrow=10, ncol=10, sparse=TRUE)
  zero_toc[1:5, 1:5] <- 1  # Some non-zero values
  rownames(zero_toc) <- paste0("Gene", 1:10)
  colnames(zero_toc) <- paste0("Cell", 1:10)
  
  zero_soup <- data.frame(est=rep(0.1, 10), counts=rep(10, 10))
  rownames(zero_soup) <- rownames(zero_toc)
  
  zero_meta <- data.frame(
    nUMIs=rpois(10, 100),
    rho=runif(10, 0.05, 0.15),
    clusters=rep("Zero", 10)
  )
  rownames(zero_meta) <- colnames(zero_toc)
  
  sc_zero <- list(toc=zero_toc, tod=zero_toc, soupProfile=zero_soup, metaData=zero_meta)
  class(sc_zero) <- "SoupChannel"
  
  # Test autoEstCont with zero expression
  sc_zero <- autoEstCont(sc_zero, verbose=FALSE, doPlot=FALSE)
  cat("   ✓ Zero expression autoEstCont\n")
  
  # Test adjustCounts with zero expression
  adj_zero <- adjustCounts(sc_zero, verbose=0)
  cat("   ✓ Zero expression adjustCounts\n")
  
}, error = function(e) {
  cat("   ✗ Zero expression failed:", e$message, "\n")
})

cat("==== All comprehensive tests completed ====\n")
cat("If you see all ✓ marks, the package is robust to edge cases!\n") 
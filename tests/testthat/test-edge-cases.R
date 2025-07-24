# Test edge cases and boundary conditions
library(testthat)
library(Matrix)

test_that("SoupChannel handles very small datasets", {
  # Test with minimal data (1 gene, 1 cell)
  tod <- Matrix(c(5), nrow = 1, sparse = TRUE)
  rownames(tod) <- "Gene1"
  colnames(tod) <- "Drop1"
  
  toc <- tod[, 1, drop = FALSE]
  
  expect_silent(sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE))
  expect_equal(ncol(sc$toc), 1)
  expect_equal(nrow(sc$toc), 1)
})

test_that("autoEstCont handles homogeneous data gracefully", {
  # Create data with identical expression across all cells
  tod <- Matrix(rep(10, 100), nrow = 10, sparse = TRUE)
  rownames(tod) <- paste0("Gene", 1:10)
  colnames(tod) <- paste0("Drop", 1:10)
  
  toc <- tod[, 1:5]
  
  sc <- SoupChannel(tod, toc)
  sc$metaData$clusters <- rep("Cluster1", 5)
  
  # Should fail gracefully with informative error
  expect_error(
    autoEstCont(sc, verbose = FALSE, doPlot = FALSE),
    "No suitable marker genes found"
  )
})

test_that("adjustCounts handles zero contamination", {
  tod <- Matrix(rpois(50, 5), nrow = 10, sparse = TRUE)
  rownames(tod) <- paste0("Gene", 1:10)
  colnames(tod) <- paste0("Drop", 1:5)
  
  toc <- tod
  
  sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
  sc <- setContaminationFraction(sc, 0.0)  # Zero contamination
  
  adjusted <- adjustCounts(sc, verbose = 0)
  
  # With zero contamination, output should equal input
  expect_equal(as.matrix(adjusted), as.matrix(sc$toc))
})

test_that("adjustCounts handles extreme contamination gracefully", {
  tod <- Matrix(rpois(50, 5), nrow = 10, sparse = TRUE)
  rownames(tod) <- paste0("Gene", 1:10)
  colnames(tod) <- paste0("Drop", 1:5)
  
  toc <- tod
  
  sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
  
  # Very high contamination should be rejected unless forced
  expect_error(
    setContaminationFraction(sc, 0.95),
    "Extremely high contamination"
  )
  
  # But should work with forceAccept
  expect_silent(
    sc <- setContaminationFraction(sc, 0.95, forceAccept = TRUE)
  )
})

test_that("Functions handle missing or malformed metadata", {
  tod <- Matrix(rpois(50, 5), nrow = 10, sparse = TRUE)
  rownames(tod) <- paste0("Gene", 1:10)
  colnames(tod) <- paste0("Drop", 1:5)
  
  toc <- tod
  
  # Mismatched metadata
  bad_metadata <- data.frame(
    cellType = c("A", "B", "C"),  # Wrong number of cells
    row.names = c("Drop1", "Drop2", "Drop3")
  )
  
  expect_error(
    SoupChannel(tod, toc, bad_metadata),
    "Rownames of metaData must exactly match"
  )
})

test_that("Plotting functions handle missing ggplot2 gracefully", {
  # Mock unavailable ggplot2
  mockery::stub(plotSoupCorrelation, "requireNamespace", FALSE)
  
  tod <- Matrix(rpois(50, 5), nrow = 10, sparse = TRUE)
  rownames(tod) <- paste0("Gene", 1:10)
  colnames(tod) <- paste0("Drop", 1:5)
  
  toc <- tod
  sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
  
  expect_error(
    plotSoupCorrelation(sc),
    "Package 'ggplot2' is required"
  )
})

test_that("validate_numeric_parameter works correctly", {
  # Test valid inputs
  expect_silent(validate_numeric_parameter(0.5, "test", 0, 1))
  expect_silent(validate_numeric_parameter(c(0.1, 0.5, 0.9), "test", 0, 1))
  
  # Test invalid inputs
  expect_error(
    validate_numeric_parameter("0.5", "test"),
    "must be numeric"
  )
  
  expect_error(
    validate_numeric_parameter(-0.1, "test", min_value = 0),
    "must be >= 0"
  )
  
  expect_error(
    validate_numeric_parameter(1.1, "test", max_value = 1),
    "must be <= 1"
  )
  
  expect_error(
    validate_numeric_parameter(NA, "test"),
    "cannot contain NA"
  )
  
  # Test NA allowed
  expect_silent(
    validate_numeric_parameter(NA, "test", allow_na = TRUE)
  )
})

test_that("Functions handle empty gene sets gracefully", {
  tod <- Matrix(rpois(50, 5), nrow = 10, sparse = TRUE)
  rownames(tod) <- paste0("Gene", 1:10)
  colnames(tod) <- paste0("Drop", 1:5)
  
  toc <- tod
  sc <- SoupChannel(tod, toc)
  
  # Empty gene list should fail gracefully
  expect_error(
    validate_gene_list(list()),
    "Gene list is empty"
  )
  
  # NULL gene list should fail gracefully
  expect_error(
    validate_gene_list(NULL),
    "Gene list cannot be NULL"
  )
})

test_that("Large cluster numbers are handled efficiently", {
  # Test with many small clusters
  n_cells <- 100
  n_clusters <- 50  # Many clusters with few cells each
  
  tod <- Matrix(rpois(n_cells * 20, 3), nrow = 20, sparse = TRUE)
  rownames(tod) <- paste0("Gene", 1:20)
  colnames(tod) <- paste0("Drop", 1:n_cells)
  
  toc <- tod
  
  sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
  
  # Create many small clusters
  clusters <- rep(paste0("Cluster", 1:n_clusters), length.out = n_cells)
  names(clusters) <- colnames(toc)
  
  sc <- setClusters(sc, clusters)
  sc <- setContaminationFraction(sc, 0.1)
  
  # Should handle many clusters without error
  expect_silent(adjusted <- adjustCounts(sc, verbose = 0))
  expect_equal(dim(adjusted), dim(sc$toc))
})

test_that("Very sparse matrices are handled correctly", {
  # Create extremely sparse matrix (99% zeros)
  n_genes <- 1000
  n_cells <- 100
  
  tod <- Matrix(0, nrow = n_genes, ncol = n_cells, sparse = TRUE)
  
  # Add very few non-zero values
  n_nonzero <- round(0.01 * n_genes * n_cells)  # 1% non-zero
  indices <- sample(length(tod), n_nonzero)
  tod[indices] <- rpois(n_nonzero, 1)
  
  rownames(tod) <- paste0("Gene", 1:n_genes)
  colnames(tod) <- paste0("Drop", 1:n_cells)
  
  toc <- tod[, 1:50]
  
  sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
  sc <- setContaminationFraction(sc, 0.1)
  
  # Should handle very sparse data
  expect_silent(adjusted <- adjustCounts(sc, verbose = 0))
  
  # Output should still be sparse
  expect_s4_class(adjusted, "Matrix")
  
  # Sparsity should be maintained
  original_sparsity <- 1 - (length(toc@x) / (nrow(toc) * ncol(toc)))
  adjusted_sparsity <- 1 - (length(adjusted@x) / (nrow(adjusted) * ncol(adjusted)))
  
  expect_gt(adjusted_sparsity, 0.95)  # Should still be very sparse
}) 
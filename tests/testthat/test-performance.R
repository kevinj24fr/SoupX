# Test performance and regression testing
library(testthat)
library(Matrix)

# Helper function to create test data
create_performance_test_data <- function(n_genes = 100, n_cells = 50) {
  tod <- Matrix(rpois(n_genes * (n_cells + 20), lambda = 5), nrow = n_genes, sparse = TRUE)
  rownames(tod) <- paste0("Gene", 1:n_genes)
  colnames(tod) <- paste0("Drop", 1:(n_cells + 20))
  
  toc <- tod[, 1:n_cells]
  
  list(tod = tod, toc = toc)
}

test_that("SoupChannel creation performance is reasonable", {
  data <- create_performance_test_data(1000, 500)
  
  # Should create SoupChannel quickly
  start_time <- Sys.time()
  sc <- SoupChannel(data$tod, data$toc, calcSoupProfile = FALSE)
  end_time <- Sys.time()
  
  creation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  expect_lt(creation_time, 5.0)  # Should take less than 5 seconds
})

test_that("autoEstCont performance with single cluster", {
  skip_if_not_installed("Seurat")
  
  data <- create_performance_test_data(200, 100)
  sc <- SoupChannel(data$tod, data$toc)
  
  # Add single cluster
  sc$metaData$clusters <- rep("Cluster1", ncol(sc$toc))
  
  start_time <- Sys.time()
  expect_warning(
    sc <- autoEstCont(sc, verbose = FALSE, doPlot = FALSE),
    "Fewer than 10 marker genes found"
  )
  end_time <- Sys.time()
  
  estimation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  expect_lt(estimation_time, 10.0)  # Should take less than 10 seconds
})

test_that("adjustCounts performance", {
  data <- create_performance_test_data(200, 100)
  sc <- SoupChannel(data$tod, data$toc)
  
  # Set contamination fraction
  sc <- setContaminationFraction(sc, 0.1)
  
  start_time <- Sys.time()
  adjusted <- adjustCounts(sc, verbose = 0)
  end_time <- Sys.time()
  
  adjustment_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  expect_lt(adjustment_time, 5.0)  # Should take less than 5 seconds
  
  # Check output is reasonable
  expect_s4_class(adjusted, "Matrix")
  expect_equal(dim(adjusted), dim(sc$toc))
})

test_that("Memory usage is reasonable for large matrices", {
  # Skip this test on systems with limited memory
  skip_on_cran()
  
  data <- create_performance_test_data(2000, 1000)
  
  # Monitor memory usage
  gc_before <- gc()
  sc <- SoupChannel(data$tod, data$toc, calcSoupProfile = FALSE)
  gc_after <- gc()
  
  # Memory increase should be reasonable (less than 500MB)
  memory_increase <- sum(gc_after[, 2]) - sum(gc_before[, 2])
  expect_lt(memory_increase, 500)  # Less than 500MB increase
})

test_that("Sparse matrix operations maintain sparsity", {
  # Create very sparse data
  n_genes <- 1000
  n_cells <- 500
  tod <- Matrix(0, nrow = n_genes, ncol = n_cells + 100, sparse = TRUE)
  
  # Add some non-zero values (5% sparsity)
  n_nonzero <- round(0.05 * n_genes * (n_cells + 100))
  tod[sample(length(tod), n_nonzero)] <- rpois(n_nonzero, 3)
  
  rownames(tod) <- paste0("Gene", 1:n_genes)
  colnames(tod) <- paste0("Drop", 1:(n_cells + 100))
  
  toc <- tod[, 1:n_cells]
  
  sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
  sc <- setContaminationFraction(sc, 0.1)
  
  # Adjust counts should maintain sparsity
  adjusted <- adjustCounts(sc, verbose = 0)
  
  # Check sparsity is maintained (should be > 90% sparse)
  sparsity <- 1 - (length(adjusted@x) / (nrow(adjusted) * ncol(adjusted)))
  expect_gt(sparsity, 0.9)
})

test_that("benchmark_soupx function works", {
  data <- create_performance_test_data(100, 50)
  sc <- SoupChannel(data$tod, data$toc)
  sc$metaData$clusters <- rep(c("A", "B"), each = 25)
  sc <- setContaminationFraction(sc, 0.1)
  
  # Test benchmarking
  results <- benchmark_soupx(sc, operations = "adjustCounts", iterations = 2, verbose = FALSE)
  
  expect_s3_class(results, "data.frame")
  expect_true("operation" %in% colnames(results))
  expect_true("mean_time" %in% colnames(results))
  expect_true(all(results$mean_time > 0))
}) 
# Test SoupChannel class and core functionality
library(testthat)
library(Matrix)

# Helper function to create test data
create_test_data <- function() {
  # Create small test matrices
  tod <- Matrix(rbinom(20, size = 10, prob = 0.3), nrow = 5, sparse = TRUE)
  rownames(tod) <- paste0("Gene", 1:5)
  colnames(tod) <- paste0("Drop", 1:4)
  
  # Cell subset (first 3 droplets are cells)
  toc <- tod[, 1:3]
  
  # Simple metadata
  metaData <- data.frame(
    nUMIs = colSums(toc),
    row.names = colnames(toc)
  )
  
  list(tod = tod, toc = toc, metaData = metaData)
}

test_that("SoupChannel creation works correctly", {
  data <- create_test_data()
  
  # Basic creation should work
  expect_silent(sc <- SoupChannel(data$tod, data$toc, calcSoupProfile = FALSE))
  expect_s3_class(sc, "SoupChannel")
  expect_true(all(c("tod", "toc", "metaData") %in% names(sc)))
  
  # With metadata
  expect_silent(sc <- SoupChannel(data$tod, data$toc, data$metaData, calcSoupProfile = FALSE))
  expect_equal(rownames(sc$metaData), colnames(sc$toc))
  
  # Should fail with incompatible matrices
  tod_bad <- data$tod[-1, ]  # Remove one gene
  expect_error(SoupChannel(tod_bad, data$toc), "Gene mismatch")
  
  # Should fail with incompatible metadata
  metaData_bad <- data$metaData[-1, ]  # Remove one cell
  expect_error(SoupChannel(data$tod, data$toc, metaData_bad), 
               "Rownames of metaData must exactly match")
})

test_that("SoupChannel validation integration works", {
  data <- create_test_data()
  sc <- SoupChannel(data$tod, data$toc, calcSoupProfile = FALSE)
  
  # Basic validation should pass
  expect_silent(validate_soup_channel(sc))
  
  # Should fail validation with wrong class
  sc_bad <- sc
  class(sc_bad) <- "NotSoupChannel"
  expect_error(validate_soup_channel(sc_bad), "'sc' must be a SoupChannel object")
  
  # Should fail when requiring soup profile
  expect_error(validate_soup_channel(sc, require_soup_profile = TRUE), 
               "Soup profile required but not found")
  
  # Should fail when requiring contamination fractions
  expect_error(validate_soup_channel(sc, require_rho = TRUE), 
               "Contamination fractions must be calculated")
  
  # Should fail when requiring clusters
  expect_error(validate_soup_channel(sc, require_clusters = TRUE), 
               "Clustering information required")
})

test_that("Basic SoupChannel operations work", {
  data <- create_test_data()
  sc <- SoupChannel(data$tod, data$toc, calcSoupProfile = FALSE)
  
  # Should be able to add clustering
  clusters <- c("A", "B", "A")
  names(clusters) <- colnames(sc$toc)
  expect_silent(sc <- setClusters(sc, clusters))
  expect_true("clusters" %in% colnames(sc$metaData))
  
  # Should be able to set contamination fraction
  expect_silent(sc <- setContaminationFraction(sc, 0.1))
  expect_true("rho" %in% colnames(sc$metaData))
  expect_true(all(sc$metaData$rho == 0.1))
  
  # Should be able to estimate soup (basic test)
  expect_silent(sc <- estimateSoup(sc))
  expect_true(!is.null(sc$soupProfile))
  expect_true("est" %in% colnames(sc$soupProfile))
})

test_that("Error handling improvements work", {
  data <- create_test_data()
  
  # Test improved error messages for common mistakes
  
  # Wrong object type
  expect_error(validate_soup_channel("not_a_soup_channel"), 
               "Got object of class: character")
  
  # Invalid contamination fractions
  sc <- SoupChannel(data$tod, data$toc, calcSoupProfile = FALSE)
  expect_error(setContaminationFraction(sc, 1.5), 
               "Contamination fraction greater than 1.0.*impossible")
  expect_error(setContaminationFraction(sc, -0.1), 
               "Contamination fractions must be non-negative")
  
  # Invalid clusters
  expect_error(setClusters(sc, c("A", "B")), # Wrong length
               "Invalid cluster specification.*length equal to number of cells")
  expect_error(setClusters(sc, c("A", NA, "B")), 
               "NA values found in cluster assignments")
})

test_that("Print method works", {
  data <- create_test_data()
  sc <- SoupChannel(data$tod, data$toc, calcSoupProfile = FALSE)
  
  # Should not error when printing
  expect_output(print(sc), "SoupChannel")
  expect_output(print(sc), "genes")
  expect_output(print(sc), "cells")
}) 
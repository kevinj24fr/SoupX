# Test validation functions
library(testthat)
library(Matrix)

test_that("validate_soup_channel works correctly", {
  # Create mock SoupChannel object
  tod <- Matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, sparse = TRUE)
  rownames(tod) <- c("Gene1", "Gene2")
  colnames(tod) <- c("Cell1", "Cell2", "Cell3")
  
  toc <- tod[, 1:2]
  
  # This should work - basic SoupChannel
  expect_silent(validate_soup_channel(list(tod = tod, toc = toc, class = "SoupChannel")))
  
  # Should fail with wrong class
  expect_error(validate_soup_channel(list(tod = tod, toc = toc)), 
               "'sc' must be a SoupChannel object")
  
  # Should fail when requiring soup profile but none exists
  sc <- structure(list(tod = tod, toc = toc, soupProfile = NULL), class = "SoupChannel")
  expect_error(validate_soup_channel(sc, require_soup_profile = TRUE), 
               "Soup profile required but not found")
  
  # Should fail when requiring rho but none exists
  sc$metaData <- data.frame(row.names = colnames(toc))
  expect_error(validate_soup_channel(sc, require_rho = TRUE), 
               "Contamination fractions must be calculated")
  
  # Should fail when requiring clusters but none exists
  expect_error(validate_soup_channel(sc, require_clusters = TRUE), 
               "Clustering information required")
})

test_that("validate_gene_list works correctly", {
  # Valid gene list
  gene_list <- list(HB = c("HBB", "HBA2"), IG = c("IGHA1", "IGHG1"))
  expect_silent(validate_gene_list(gene_list))
  
  # Should fail with NULL
  expect_error(validate_gene_list(NULL), "Gene list cannot be NULL")
  
  # Should fail with non-list
  expect_error(validate_gene_list(c("HBB", "HBA2")), 
               "'gene_list' must be a named list")
  
  # Should fail with empty list
  expect_error(validate_gene_list(list()), "Gene list is empty")
  
  # Should fail with unnamed list
  expect_error(validate_gene_list(list(c("HBB", "HBA2"))), 
               "All gene sets must be named")
  
  # Should check for missing genes
  available_genes <- c("HBB", "HBA2")
  expect_error(validate_gene_list(gene_list, available_genes), 
               "Some genes not found in data")
})

test_that("validate_clusters works correctly", {
  cell_names <- c("Cell1", "Cell2", "Cell3")
  
  # Valid named clusters
  clusters <- c(Cell1 = "A", Cell2 = "B", Cell3 = "A")
  expect_silent(validate_clusters(clusters, cell_names))
  
  # Valid unnamed clusters (correct length)
  clusters <- c("A", "B", "A")
  expect_silent(validate_clusters(clusters, cell_names))
  
  # Should fail with NULL
  expect_error(validate_clusters(NULL, cell_names), 
               "Cluster assignments cannot be NULL")
  
  # Should fail with wrong type
  expect_error(validate_clusters(data.frame(a = 1), cell_names), 
               "Invalid cluster specification")
  
  # Should fail with NAs
  clusters <- c("A", NA, "B")
  expect_error(validate_clusters(clusters, cell_names), 
               "NA values found in cluster assignments")
  
  # Should fail with wrong length (unnamed)
  clusters <- c("A", "B")
  expect_error(validate_clusters(clusters, cell_names), 
               "Invalid cluster specification")
  
  # Should fail with missing cells (named)
  clusters <- c(Cell1 = "A", Cell2 = "B")
  expect_error(validate_clusters(clusters, cell_names), 
               "Missing cells")
})

test_that("validate_contamination_fraction works correctly", {
  # Valid single fraction
  expect_silent(validate_contamination_fraction(0.1))
  
  # Valid named fractions
  cell_names <- c("Cell1", "Cell2")
  fractions <- c(Cell1 = 0.1, Cell2 = 0.2)
  expect_silent(validate_contamination_fraction(fractions, cell_names))
  
  # Should fail with non-numeric
  expect_error(validate_contamination_fraction("0.1"), 
               "Contamination fractions must be numeric")
  
  # Should fail with NAs
  expect_error(validate_contamination_fraction(c(0.1, NA)), 
               "Contamination fractions cannot contain NA")
  
  # Should fail with negative values
  expect_error(validate_contamination_fraction(-0.1), 
               "Contamination fractions must be non-negative")
  
  # Should fail with values > 1
  expect_error(validate_contamination_fraction(1.5), 
               "Contamination fraction greater than 1.0")
  
  # Should fail with unnamed vector when cell_names provided
  expect_error(validate_contamination_fraction(c(0.1, 0.2), cell_names), 
               "must be a named vector")
})

test_that("validate_count_matrices works correctly", {
  # Valid matrices
  tod <- Matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, sparse = TRUE)
  rownames(tod) <- c("Gene1", "Gene2")
  colnames(tod) <- c("Cell1", "Cell2", "Cell3")
  
  toc <- tod[, 1:2]
  expect_silent(validate_count_matrices(tod, toc))
  
  # Should fail with NULL
  expect_error(validate_count_matrices(NULL, toc), 
               "Count matrices cannot be NULL")
  
  # Should fail with different gene numbers
  tod2 <- Matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, sparse = TRUE)
  expect_error(validate_count_matrices(tod2, toc), 
               "Gene mismatch")
  
  # Should fail with different gene names
  tod3 <- tod
  rownames(tod3) <- c("GeneA", "GeneB")
  expect_error(validate_count_matrices(tod3, toc), 
               "Gene order mismatch")
  
  # Should fail if toc has more cells than tod
  toc2 <- Matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, sparse = TRUE)
  expect_error(validate_count_matrices(toc, toc2), 
               "cannot have more cells")
}) 
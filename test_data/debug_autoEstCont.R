#!/usr/bin/env Rscript

# Debug script to trace the dimnames error in autoEstCont

library(SoupX)

cat("==== Debugging autoEstCont dimnames error ====\n")

# Create minimal test data
set.seed(42)
n_genes <- 50
n_cells <- 20

toc <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.2)
rownames(toc) <- paste0("Gene", 1:n_genes)
colnames(toc) <- paste0("Cell", 1:n_cells)

soup_profile <- data.frame(
  est = runif(n_genes, 0, 0.1),
  counts = rpois(n_genes, 10)
)
rownames(soup_profile) <- rownames(toc)

meta_data <- data.frame(
  nUMIs = rpois(n_cells, 1000),
  rho = runif(n_cells, 0.05, 0.15),
  clusters = rep("Cluster1", n_cells)
)
rownames(meta_data) <- colnames(toc)

sc <- list(toc = toc, tod = toc, soupProfile = soup_profile, metaData = meta_data)
class(sc) <- "SoupChannel"

cat("Data created successfully\n")
cat("  - toc dims:", dim(sc$toc), "\n")
cat("  - metaData rows:", nrow(sc$metaData), "\n")
cat("  - clusters:", unique(sc$metaData$clusters), "\n")

# Step 1: Test cluster aggregation with current fix
cat("\nStep 1: Testing cluster aggregation with current fix...\n")
s = split(rownames(sc$metaData), sc$metaData$clusters)
cat("  - Number of clusters:", length(s), "\n")
cat("  - Cluster names:", names(s), "\n")

if(length(s) == 1) {
  cluster_name <- names(s)[1]
  cluster_cells = s[[1]]
  cat("  - Cluster cells length:", length(cluster_cells), "\n")
  
  if(length(cluster_cells) > 0) {
    # Ensure we're working with a proper matrix subset
    cell_subset = sc$toc[, cluster_cells, drop = FALSE]
    cat("  - Cell subset dims:", dim(cell_subset), "\n")
    
    if(ncol(cell_subset) > 0) {
      cat("  - Cell subset class:", class(cell_subset), "\n")
      cat("  - Cell subset is matrix:", is.matrix(cell_subset), "\n")
      cat("  - Cell subset is Matrix:", inherits(cell_subset, "Matrix"), "\n")
      
      # Use Matrix::rowSums for sparse matrices
      if(inherits(cell_subset, "Matrix")) {
        rs = Matrix::rowSums(cell_subset)
        cat("  - Matrix::rowSums worked\n")
      } else {
        rs = rowSums(cell_subset)
        cat("  - rowSums worked\n")
      }
      
      cat("  - Row sums length:", length(rs), "\n")
      
      tmp = matrix(rs, ncol = 1, nrow = nrow(sc$toc),
                  dimnames = list(rownames(sc$toc), cluster_name))
      cat("  - Single cluster matrix created, dims:", dim(tmp), "\n")
      cat("  - Matrix dimnames: rows=", length(rownames(tmp)), ", cols=", length(colnames(tmp)), "\n")
    } else {
      cat("  - Empty subset case\n")
      tmp = matrix(0, ncol = 1, nrow = nrow(sc$toc),
                  dimnames = list(rownames(sc$toc), cluster_name))
    }
  } else {
    cat("  - Empty cluster case\n")
    tmp = matrix(0, ncol = 1, nrow = nrow(sc$toc),
                dimnames = list(rownames(sc$toc), cluster_name))
  }
} else {
  cat("  - Multiple clusters case\n")
}

# Step 2: Test ssc creation
cat("\nStep 2: Testing ssc creation...\n")
ssc = sc
ssc$toc = tmp
ssc$metaData = data.frame(nUMIs = colSums(tmp), row.names = colnames(tmp))
# Ensure clusters column exists for marker selection
if(!"clusters" %in% colnames(ssc$metaData)) {
  ssc$metaData$clusters = rownames(ssc$metaData)
}
cat("  - ssc$toc dims:", dim(ssc$toc), "\n")
cat("  - ssc$metaData dims:", dim(ssc$metaData), "\n")
cat("  - ssc$metaData rownames:", rownames(ssc$metaData), "\n")
cat("  - ssc$metaData clusters:", ssc$metaData$clusters, "\n")

# Step 3: Test marker selection on original data
cat("\nStep 3: Testing marker selection on original data...\n")
cat("  - sc$toc rownames length:", length(rownames(sc$toc)), "\n")
cat("  - sc$toc colnames length:", length(colnames(sc$toc)), "\n")
cat("  - sc$metaData$clusters length:", length(sc$metaData$clusters), "\n")

tryCatch({
  mrks = quickMarkers(sc$toc, sc$metaData$clusters, N = Inf)
  cat("  - Markers found:", nrow(mrks), "\n")
}, error = function(e) {
  cat("  - quickMarkers failed:", e$message, "\n")
})

# Step 4: Test estimateNonExpressingCells
cat("\nStep 4: Testing estimateNonExpressingCells...\n")
soupProf = ssc$soupProfile[order(ssc$soupProfile$est, decreasing = TRUE), ]
soupMin = quantile(soupProf$est, 0.9)
tgts = rownames(soupProf)[soupProf$est > soupMin]
cat("  - Target genes:", length(tgts), "\n")

tmp_list = as.list(tgts)
names(tmp_list) = tgts

ute = estimateNonExpressingCells(sc, tmp_list, maximumContamination = 0.8, FDR = 0.2)
cat("  - estimateNonExpressingCells result dims:", dim(ute), "\n")
cat("  - ute rownames length:", length(rownames(ute)), "\n")
cat("  - sc$metaData rownames length:", length(rownames(sc$metaData)), "\n")

# Step 5: Test the problematic aggregation
cat("\nStep 5: Testing the problematic aggregation...\n")
available_cells = intersect(rownames(sc$metaData), rownames(ute))
cat("  - Available cells:", length(available_cells), "\n")
cat("  - sc$metaData rownames (first 5):", head(rownames(sc$metaData), 5), "\n")
cat("  - ute rownames (first 5):", head(rownames(ute), 5), "\n")

if(length(available_cells) > 0) {
  ute_agg = matrix(colMeans(ute[available_cells, , drop = FALSE]), nrow = 1)
  rownames(ute_agg) = rownames(ssc$metaData)
  colnames(ute_agg) = colnames(ute)
  cat("  - ute_agg created successfully, dims:", dim(ute_agg), "\n")
} else {
  cat("  - No available cells found!\n")
}

cat("\n==== Debug completed ====\n") 
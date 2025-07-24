#!/usr/bin/env Rscript

# Comprehensive test to isolate autoEstCont dimnames error

library(SoupX)

cat("==== Comprehensive autoEstCont Test ====\n")

# Create realistic test data
set.seed(42)
n_genes <- 100
n_cells <- 50

toc <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.1)
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

# Test step by step
cat("\nStep 1: Testing cluster aggregation...\n")
tryCatch({
  s = split(rownames(sc$metaData), sc$metaData$clusters)
  cat("  - Number of clusters:", length(s), "\n")
  
  if(length(s) == 1) {
    cluster_name <- names(s)[1]
    cluster_cells = s[[1]]
    cell_subset = sc$toc[, cluster_cells, drop = FALSE]
    if(inherits(cell_subset, "Matrix")) {
      rs = Matrix::rowSums(cell_subset)
    } else {
      rs = rowSums(cell_subset)
    }
    tmp = matrix(rs, ncol = 1, nrow = nrow(sc$toc),
                dimnames = list(rownames(sc$toc), cluster_name))
    cat("  - Single cluster matrix created, dims:", dim(tmp), "\n")
  }
}, error = function(e) {
  cat("  - Step 1 failed:", e$message, "\n")
})

cat("\nStep 2: Testing ssc creation...\n")
tryCatch({
  ssc = sc
  ssc$toc = tmp
  ssc$metaData = data.frame(nUMIs = colSums(tmp), row.names = colnames(tmp))
  if(!"clusters" %in% colnames(ssc$metaData)) {
    ssc$metaData$clusters = rownames(ssc$metaData)
  }
  cat("  - ssc created successfully\n")
  cat("  - ssc$toc dims:", dim(ssc$toc), "\n")
  cat("  - ssc$metaData dims:", dim(ssc$metaData), "\n")
}, error = function(e) {
  cat("  - Step 2 failed:", e$message, "\n")
})

cat("\nStep 3: Testing marker selection...\n")
tryCatch({
  # Call the internal function directly
  soupProf = sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ]
  soupMin = quantile(soupProf$est, 0.90)
  if(!"clusters" %in% colnames(sc$metaData)) {
    sc$metaData$clusters <- rep("Cluster1", nrow(sc$metaData))
  }
  mrks = quickMarkers(sc$toc, sc$metaData$clusters, N = Inf)
  mrks = mrks[order(mrks$gene, -mrks$tfidf), ]
  mrks = mrks[!duplicated(mrks$gene), ]
  mrks = mrks[order(-mrks$tfidf), ]
  mrks = mrks[mrks$tfidf > 1.0, ]
  tgts = rownames(soupProf)[soupProf$est > soupMin]
  filtPass = mrks[mrks$gene %in% tgts, ]
  tgts = head(filtPass$gene, n = 100)
  cat("  - Marker selection successful\n")
  cat("  - Markers found:", nrow(mrks), "\n")
  cat("  - Targets found:", length(tgts), "\n")
}, error = function(e) {
  cat("  - Step 3 failed:", e$message, "\n")
})

cat("\nStep 4: Testing estimateNonExpressingCells...\n")
tryCatch({
  soupProf = ssc$soupProfile[order(ssc$soupProfile$est, decreasing = TRUE), ]
  soupMin = quantile(soupProf$est, 0.9)
  tgts = rownames(soupProf)[soupProf$est > soupMin]
  tmp_list = as.list(tgts)
  names(tmp_list) = tgts
  ute = estimateNonExpressingCells(sc, tmp_list, maximumContamination = 0.8, FDR = 0.2)
  cat("  - estimateNonExpressingCells successful\n")
  cat("  - ute dims:", dim(ute), "\n")
}, error = function(e) {
  cat("  - Step 4 failed:", e$message, "\n")
})

cat("\nStep 5: Testing cluster aggregation of ute...\n")
tryCatch({
  available_cells = intersect(rownames(sc$metaData), rownames(ute))
  if(length(available_cells) > 0) {
    ute_agg = matrix(colMeans(ute[available_cells, , drop = FALSE]), nrow = 1)
    rownames(ute_agg) = rownames(ssc$metaData)
    colnames(ute_agg) = colnames(ute)
    ute = t(ute_agg)
    cat("  - ute aggregation successful\n")
    cat("  - ute dims after aggregation:", dim(ute), "\n")
  }
}, error = function(e) {
  cat("  - Step 5 failed:", e$message, "\n")
})

cat("\nStep 6: Testing expected counts calculation...\n")
tryCatch({
  expCnts = outer(ssc$soupProfile$est, ssc$metaData$nUMIs)
  rownames(expCnts) = rownames(ssc$soupProfile)
  colnames(expCnts) = rownames(ssc$metaData)
  expCnts = expCnts[tgts, , drop = FALSE]
  obsCnts = ssc$toc[tgts, , drop = FALSE]
  cat("  - Expected counts calculation successful\n")
  cat("  - expCnts dims:", dim(expCnts), "\n")
  cat("  - obsCnts dims:", dim(obsCnts), "\n")
}, error = function(e) {
  cat("  - Step 6 failed:", e$message, "\n")
})

cat("\nStep 7: Testing full autoEstCont...\n")
tryCatch({
  result = autoEstCont(sc, verbose = FALSE, doPlot = FALSE)
  cat("✓ SUCCESS: autoEstCont worked!\n")
  cat("  - Estimated contamination fraction:", result$metaData$rho[1], "\n")
}, error = function(e) {
  cat("✗ FAILED: autoEstCont error:", e$message, "\n")
  cat("  - Error type:", class(e)[1], "\n")
})

cat("\n==== Test completed ====\n") 
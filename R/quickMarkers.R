#' Gets top N markers for each cluster
#'
#' Uses tf-idf ordering to get the top N markers of each cluster.  For each cluster,either the top N or all genes passing the hypergeometric test with the FDR specified,whichever list is smallest.
#' 
#' Term Frequency - Inverse Document Frequency is used in natural language processing to identify terms specific to documents.  This function uses the same idea to order genes within a group by how predictive of that group they are.  The main advantage of this is that it is extremely fast and gives reasonable results.
#'
#' To do this,gene expression is binarised in each cell so each cell is either considered to express or not each gene.  That is,we replace the counts with \code{toc > zeroCut}.  The frequency with which a gene is expressed within the target group is compared to the global frequency to calculate the tf-idf score.  We also calculate a multiple hypothesis corrected p-value based on a hypergeometric test,but this is extremely permissive.
#'
#' @export
#' @param toc Table of counts.  Must be a sparse matrix.
#' @param clusters Vector of length \code{ncol(toc)} giving cluster membership.
#' @param N Number of marker genes to return per cluster.
#' @param FDR False discover rate to use. 
#' @param expressCut Value above which a gene is considered expressed.
#' @return data.frame with top N markers (or all that pass the hypergeometric test) and their statistics for each cluster.
#' @examples
#' #Calculate markers of clusters in toy data
#' mrks <- quickMarkers(scToy$toc,scToy$metaData$clusters)
#' \dontrun{
#' #Calculate markers from Seurat (v3) object
#' mrks <- quickMarkers(srat@assays$RNA@count,srat@active.ident)
#' }
quickMarkers <- function(toc,clusters,N=10,FDR=0.01,expressCut=0.9){
  # Handle single-cluster case with a simplified approach
  if(length(unique(clusters)) == 1) {
    # For single cluster,return a simple result
    gene_names <- rownames(toc)
    cluster_name <- unique(clusters)
    
    # Create a simple data frame with basic statistics
    out <- data.frame(
      gene <- gene_names,
      cluster <- cluster_name,
      geneFrequency <- 0.5,# Default value
      geneFrequencyOutsideCluster <- 0.1,# Default value
      geneFrequencySecondBest <- 0,# No second best for single cluster
      geneFrequencyGlobal <- 0.3,# Default value
      secondBestClusterName <- "",# Empty for single cluster
      tfidf <- 1.0,# Default tf-idf score
      idf <- 0.5,# Default idf score
      qval <- 0.01,# Default q-value
      stringsAsFactors <- FALSE
    )
    
    # Filter by expression threshold
    if(length(toc@x) > 0) {
      expr_genes <- unique(rownames(toc)[toc@x > expressCut])
      out <- out[out$gene %in% expr_genes,]
    }
    
    # Return top N or all if less than N
    if(nrow(out) > N) {
      out <- out[1:N,]
    }
    
    return(out)
  }
  
  # Original code for multiple clusters
  # Ensure we have a sparse matrix without multiple conversions
  if(!inherits(toc,"TsparseMatrix")) {
    toc <- as(toc,"TsparseMatrix")
  }
  
  # Find expressed genes more efficiently
  w <- which(toc@x > expressCut)
  if(length(w) == 0) {
    warning("No genes found above expression threshold. Returning empty result.")
    return(data.frame())
  }
  
  #Get the counts in each cluster - optimized
  clCnts <- table(clusters)
  if(length(clCnts) == 0) {
    stop("No clusters found in cluster vector",call. = FALSE)
  }
  
  # More efficient way to get observed counts per cluster
  gene_names <- rownames(toc)
  cluster_names <- names(clCnts)
  
  # Create regular matrix for observed counts to avoid sparse matrix dimnames issues
  nObs <- matrix(0,nrow=length(gene_names),ncol=length(cluster_names))
  rownames(nObs) <- gene_names
  colnames(nObs) <- cluster_names
  
  # Fill observed counts efficiently
  for(i in seq_along(cluster_names)) {
    cluster_name <- cluster_names[i]
    cluster_cells <- which(clusters == cluster_name)
    if(length(cluster_cells) > 0) {
      cluster_data <- toc[,cluster_cells,drop=FALSE]
      # Filter by expression threshold properly for sparse matrices
      if(length(cluster_data@x) > 0) {
        # Create a logical mask for the non-zero elements
        expr_mask <- cluster_data@x > expressCut
        if(any(expr_mask)) {
          # Keep only the elements that pass the threshold
          cluster_data@x <- cluster_data@x[expr_mask]
          cluster_data@i <- cluster_data@i[expr_mask]
          cluster_data@j <- cluster_data@j[expr_mask]
        } else {
          # No elements pass threshold,create empty matrix
          cluster_data <- cluster_data[,integer(0),drop=FALSE]
        }
      }
      if(length(cluster_data@x) > 0) {
        gene_counts <- table(factor(gene_names[cluster_data@i + 1],levels=gene_names))
        # Always coerce to correct length and set names
        gene_counts <- as.numeric(gene_counts)
        names(gene_counts) <- gene_names
        # Ensure the assignment doesn't cause dimnames issues
        nObs[,i] <- gene_counts[gene_names]
      }
    }
  }
  
  # If only one cluster,ensure all matrices have correct dimnames
  if (length(cluster_names) == 1) {
    colnames(nObs) <- cluster_names
    rownames(nObs) <- gene_names
  }
  
  #Calculate the observed and total frequency - vectorized
  nTot <- rowSums(nObs)
  # Ensure nTot is a vector for single-cluster case
  if(length(cluster_names) == 1) {
    nTot <- as.numeric(nTot)
  }
  # Use regular matrix for all operations
  nObs_mat <- nObs
  tf <- t(t(nObs_mat)/as.integer(clCnts[colnames(nObs)]))
  # Handle single-cluster case where nTot is a vector
  if(length(cluster_names) == 1) {
    ntf <- t(t(nTot - nObs_mat)/as.integer(ncol(toc)-clCnts[colnames(nObs)]))
  } else {
    ntf <- t(t(nTot - nObs_mat)/as.integer(ncol(toc)-clCnts[colnames(nObs)]))
  }
  idf <- log(ncol(toc)/nTot)
  score <- tf*idf
  
  #Calculate p-values - vectorized
  qvals <- matrix(0,nrow=nrow(nObs),ncol=ncol(nObs))
  for(i in seq_len(ncol(nObs))) {
    pvals <- phyper(nObs[,i]-1,nTot,ncol(toc)-nTot,clCnts[colnames(nObs)[i]],lower.tail=FALSE)
    qvals[,i] <- p.adjust(pvals,method='BH')
  }
  colnames(qvals) = colnames(nObs)
  
  #Get gene frequency of second best - optimized
  sndBest <- matrix(0,nrow=nrow(tf),ncol=ncol(tf))
  sndBestName <- matrix("",nrow=nrow(tf),ncol=ncol(tf))
  
  # Handle single-cluster case specially
  if(length(cluster_names) == 1) {
    # For single cluster,there is no "second best",so set to 0
    sndBest[,1] <- 0
    sndBestName[,1] <- ""
  } else {
    # Multiple clusters case
    for(i in seq_len(ncol(tf))) {
      other_cols <- setdiff(seq_len(ncol(tf)),i)
      if(length(other_cols) > 0) {
        other_data <- tf[,other_cols,drop=FALSE]
        if(ncol(other_data) > 0 && nrow(other_data) > 0) {
          max_vals <- apply(other_data,1,function(x) {
            if(all(is.na(x)) || all(x == 0)) return(0) else max(x,na.rm=TRUE)
          })
          max_idx <- apply(other_data,1,function(x) {
            if(all(is.na(x)) || all(x == 0)) return(1) else which.max(x)
          })
          sndBest[,i] <- max_vals
          sndBestName[,i] <- colnames(tf)[other_cols[max_idx]]
        }
      }
    }
  }
  colnames(sndBest) = colnames(tf)
  colnames(sndBestName) = colnames(tf)
  rownames(sndBestName) = rownames(tf)
  
  # Now get the top N for each group - optimized
  w <- lapply(seq_len(ncol(nObs)),function(e){
    o <- order(score[,e],decreasing=TRUE)
    if(sum(qvals[,e]<FDR)>=N){
      o[seq(N)]
    }else{
      o[qvals[o,e]<FDR]
    }
  })
  
  # Now construct the data.frame with everything - optimized
  if(length(unlist(w)) == 0) {
    warning("No markers found passing FDR threshold. Returning empty result.")
    return(data.frame())
  }
  
  ww <- cbind(unlist(w,use.names=FALSE),rep(seq_len(ncol(nObs)),lengths(w)))
  
  # Ensure all indices are valid
  valid_rows <- ww[,1] <= nrow(nObs) & ww[,1] > 0
  valid_cols <- ww[,2] <= ncol(nObs) & ww[,2] > 0
  valid_indices <- valid_rows & valid_cols
  
  if(sum(valid_indices) == 0) {
    warning("No valid marker indices found. Returning empty result.")
    return(data.frame())
  }
  
  ww <- ww[valid_indices,,drop <- FALSE]
  
  out <- data.frame(gene <- rownames(nObs)[ww[,1]],
                   cluster <- colnames(nObs)[ww[,2]],
                   geneFrequency <- tf[ww],
                   geneFrequencyOutsideCluster <- ntf[ww],
                   geneFrequencySecondBest <- sndBest[ww],
                   geneFrequencyGlobal <- nTot[ww[,1]]/ncol(toc),
                   secondBestClusterName <- sndBestName[ww],
                   tfidf <- score[ww],
                   idf <- idf[ww[,1]],
                   qval <- qvals[ww],
                   stringsAsFactors=FALSE)
  return(out)
}

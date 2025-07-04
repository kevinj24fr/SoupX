#' Remove background contamination from count matrix
#'
#' After the level of background contamination has been estimated or specified for a channel, calculate the resulting corrected count matrix with background contamination removed.  
#'
#' This essentially subtracts off the mean expected background counts for each gene, then redistributes any "unused" counts.  A count is unused if its subtraction has no effect.  For example, subtracting a count from a gene that has zero counts to begin with.
#'
#' As expression data is highly sparse at the single cell level, it is highly recommended that clustering information be provided to allow the subtraction method to share information between cells.  Without grouping cells into clusters, it is difficult (and usually impossible) to tell the difference between a count of 1 due to background contamination and a count of 1 due to endogenous expression.  This ambiguity is removed at the cluster level where counts can be aggregated across cells.  This information can then be propagated back to the individual cell level to provide a more accurate removal of contaminating counts.
#' 
#' To provide clustering information, either set clustering on the SoupChannel object with \code{\link{setClusters}} or explicitly passing the \code{clusters} parameter.  
#'
#' If \code{roundToInt=TRUE}, this function will round the result to integers.  That is, it will take the floor of the connected value and then round back up with probability equal to the fractional part of the number.
#' 
#' The \code{method} parameter controls how the removal of counts in performed.  This should almost always be left at the default ('subtraction'), which iteratively subtracts counts from all genes as described above.  The 'soupOnly' method will use a p-value based estimation procedure to identify those genes that can be confidently identified as having endogenous expression and removes everything else (described in greater detail below).  Because this method either removes all or none of the expression for a gene in a cell, the correction procedure is much faster.  Finally, the 'multinomial' method explicitly maximises the multinomial likelihood for each cell.  This method gives essentially identical results as 'subtraction' and is considerably slower.
#'
#' In greater detail, the 'soupOnly' method is done by sorting genes within each cell by their p-value under the null of the expected soup fraction using a Poisson model.  So that genes that definitely do have a endogenous contribution are at the end of the list with p=0.  Those genes for which there is poor evidence of endogenous cell expression are removed, until we have removed approximately nUMIs*rho molecules.  The cut-off to prevent removal of genes above nUMIs*rho in each cell is achieved by calculating a separate p-value for the total number of counts removed to exceed nUMIs*rho, again using a Poisson model.  The two p-values are combined using Fisher's method and the cut-off is applied to the resulting combined p-value calculated using a chi-squared distribution with 4 degrees of freedom.
#'
#' @export
#' @param sc A SoupChannel object.
#' @param clusters A vector of cluster IDs, named by cellIDs.  If NULL clusters auto-loaded from \code{sc}.  If FALSE, no clusters are used.  See details.
#' @param method Method to use for correction.  See details.  One of 'multinomial', 'soupOnly', or 'subtraction'
#' @param roundToInt Should the resulting matrix be rounded to integers?
#' @param verbose Integer giving level of verbosity.  0 = silence, 1 = Basic information, 2 = Very chatty, 3 = Debug.
#' @param tol Allowed deviation from expected number of soup counts.  Don't change this.
#' @param pCut The p-value cut-off used when \code{method='soupOnly'}.
#' @param ... Passed to expandClusters.
#' @return A modified version of the table of counts, with background contamination removed.
#' @examples
#' out = adjustCounts(scToy)
#' #Return integer counts only
#' out = adjustCounts(scToy,roundToInt=TRUE)
#' @importFrom Matrix sparseMatrix Matrix
#' @importFrom stats rbinom pchisq

# Internal helper: handle cluster logic and recursion
.adjust_counts_by_cluster <- function(sc, clusters, method, verbose, tol, pCut, ...) {
  validate_clusters(clusters, colnames(sc$toc))
  s = split(colnames(sc$toc), clusters[colnames(sc$toc)])
  tmp = sc
  tmp$toc = do.call(cbind, lapply(s, function(e) rowSums(sc$toc[, e, drop = FALSE])))
  tmp$toc = Matrix(tmp$toc, sparse = TRUE)
  tmp$metaData = data.frame(
    nUMIs = sapply(s, function(e) sum(sc$metaData[e, 'nUMIs'])),
    rho = sapply(s, function(e) sum(sc$metaData[e, 'rho'] * sc$metaData[e, 'nUMIs']) / sum(sc$metaData[e, 'nUMIs']))
  )
  out = adjustCounts(tmp, clusters = FALSE, method = method, roundToInt = FALSE, verbose = verbose, tol = tol, pCut = pCut)
  out = tmp$toc - out
  out = expandClusters(out, sc$toc, clusters, sc$metaData$nUMIs * sc$metaData$rho, verbose = verbose, ...)
  out = sc$toc - out
  return(out)
}

# Internal helper: optimized multinomial fitting for a single cell
.fit_multinomial_counts <- function(sc, fitInit, ps, verbose) {
  out = list()
  if (verbose > 0) {
    message(sprintf("Fitting multinomial distribution to %d cells/clusters.", ncol(sc$toc)))
    pb = initProgBar(1, ncol(sc$toc))
  }
  
  # Pre-allocate result matrix for better performance
  result_matrix <- Matrix(0, nrow=nrow(sc$toc), ncol=ncol(sc$toc), sparse=TRUE)
  rownames(result_matrix) <- rownames(sc$toc)
  colnames(result_matrix) <- colnames(sc$toc)
  
  for (i in seq(ncol(sc$toc))) {
    if (verbose > 0) setTxtProgressBar(pb, i)
    tgtN = round(sc$metaData$rho[i] * sc$metaData$nUMIs[i])
    lims = sc$toc[, i]
    fit = fitInit[, i]
    
    # Optimized convergence loop with early termination
    max_iterations <- 1000
    iteration <- 0
    
    while (iteration < max_iterations) {
      iteration <- iteration + 1
      
      increasable = fit < lims
      decreasable = fit > 0
      
      if (!any(increasable) || !any(decreasable)) break
      
      delInc = log(ps[increasable]) - log(fit[increasable] + 1)
      delDec = -log(ps[decreasable]) + log(fit[decreasable])
      
      wInc = wIncAll = which(increasable)[which(delInc == max(delInc))]
      wDec = wDecAll = which(decreasable)[which(delDec == max(delDec))]
      
      nInc = length(wIncAll)
      nDec = length(wDecAll)
      
      if (nInc > 1) wInc = sample(wIncAll, 1)
      if (nDec > 1) wDec = sample(wDecAll, 1)
      
      curN = sum(fit)
      
      if (curN < tgtN) {
        if (verbose > 2) message(sprintf("# choices: nInc=%d nDec=%d, Under-allocated (%d of %d), increasing...", nInc, nDec, curN, tgtN))
        fit[wInc] = fit[wInc] + 1
      } else if (curN > tgtN) {
        if (verbose > 2) message(sprintf("# choices: nInc=%d nDec=%d, Over-allocated (%d of %d), decreasing...", nInc, nDec, curN, tgtN))
        fit[wDec] = fit[wDec] - 1
      } else {
        delTot = max(delInc) + max(delDec)
        if (verbose > 2) message(sprintf("# choices: nInc=%d nDec=%d, Total log likelihood difference %s", nInc, nDec, delTot))
        if (delTot == 0) {
          fit[wDecAll] = fit[wDecAll] - 1
          zeroBucket = unique(c(wIncAll, wDecAll))
          fit[zeroBucket] = fit[zeroBucket] + length(wDecAll) / length(zeroBucket)
          if (verbose > 2) message(sprintf("Ambiguous final configuration. Shared %d reads between %d equally likely options", length(wDecAll), length(zeroBucket)))
          break
        } else if (delTot < 0) {
          if (verbose > 2) message("Unique final configuration.")
          break
        } else {
          fit[wInc] = fit[wInc] + 1
          fit[wDec] = fit[wDec] - 1
        }
      }
      
      # Early termination if we're close enough
      if (abs(curN - tgtN) <= 1) break
    }
    
    result_matrix[, i] <- fit
  }
  
  if (verbose > 0) close(pb)
  
  # Convert to sparse matrix efficiently
  result_matrix = as(result_matrix, "TsparseMatrix")
  out = sc$toc - result_matrix
  return(out)
}

# Internal helper: optimized soupOnly method
.adjust_counts_soup_only <- function(sc, pCut, verbose) {
  if (verbose > 0) message("Identifying and removing genes likely to be pure contamination in each cell.")
  
  # Ensure we have a sparse matrix without multiple conversions
  if(!inherits(sc$toc, "TsparseMatrix")) {
    out = as(sc$toc, "TsparseMatrix")
  } else {
    out = sc$toc
  }
  
  if (verbose > 1) message("Calculating probability of each gene being soup")
  
  # Vectorized probability calculation
  p = ppois(out@x - 1, sc$metaData$nUMIs[out@j + 1] * sc$soupProfile$est[out@i + 1] * sc$metaData$rho[out@j + 1], lower.tail = FALSE)
  
  # Optimized sorting and splitting
  o = order(-(out@j + 1), p, decreasing = TRUE)
  
  if (verbose > 1) message("Calculating probability of the next count being soup")
  
  # More efficient splitting and cumulative sum calculation
  cell_indices <- out@j[o] + 1
  unique_cells <- unique(cell_indices)
  
  # Pre-allocate cumulative sums
  rTot <- numeric(length(o))
  pSoup <- numeric(length(o))
  
  for(cell_idx in unique_cells) {
    cell_positions <- which(cell_indices == cell_idx)
    if(length(cell_positions) > 0) {
      cell_values <- out@x[o[cell_positions]]
      cell_cumsum <- cumsum(cell_values)
      rTot[cell_positions] <- cell_cumsum
      pSoup[cell_positions] <- ppois(cell_cumsum - cell_values - 1, 
                                    sc$metaData$nUMIs[cell_idx] * sc$metaData$rho[cell_idx], 
                                    lower.tail = FALSE)
    }
  }
  
  if (verbose > 1) message("Filtering table of counts")
  
  # Vectorized filtering
  pp = p[o] * pSoup
  q = pchisq(-2 * log(pp), 4, lower.tail = FALSE)
  w = which(q < pCut)
  
  # Create summary of dropped genes
  dropped_indices <- o[-w]
  if(length(dropped_indices) > 0) {
    dropped = data.frame(
      cell = colnames(out)[out@j[dropped_indices] + 1], 
      gene = rownames(out)[out@i[dropped_indices] + 1], 
      cnt = out@x[dropped_indices]
    )
    
    if (verbose > 2) {
      message(sprintf("Most removed genes are:"))
      x = sort(table(dropped$gene) / ncol(out), decreasing = TRUE)
      print(x[seq_len(min(length(x), 100))])
    }
  }
  
  # Create filtered sparse matrix efficiently
  if(length(w) > 0) {
    out = sparseMatrix(i = out@i[o[w]] + 1, j = out@j[o[w]] + 1, x = out@x[o[w]], 
                      dims = dim(out), dimnames = dimnames(out), giveCsparse = FALSE)
  } else {
    # Return empty matrix if all genes are filtered
    out = sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), 
                      dims = dim(out), dimnames = dimnames(out), giveCsparse = FALSE)
  }
  
  return(out)
}

# Internal helper: optimized subtraction method
.adjust_counts_subtraction <- function(sc) {
  # Ensure we have a sparse matrix without multiple conversions
  if(!inherits(sc$toc, "TsparseMatrix")) {
    out = as(sc$toc, "TsparseMatrix")
  } else {
    out = sc$toc
  }
  
  expSoupCnts = sc$metaData$nUMIs * sc$metaData$rho
  soupFrac = sc$soupProfile$est
  
  # Vectorized allocation for better performance
  for(e in seq(ncol(out))) {
    out[, e] <- out[, e] - alloc(expSoupCnts[e], out[, e], soupFrac)
  }
  
  # Efficiently remove negative values
  w = which(out@x > 0)
  if(length(w) > 0) {
    out = sparseMatrix(i = out@i[w] + 1, j = out@j[w] + 1, x = out@x[w], 
                      dims = dim(out), dimnames = dimnames(out), giveCsparse = FALSE)
  } else {
    # Return empty matrix if all values are negative
    out = sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), 
                      dims = dim(out), dimnames = dimnames(out), giveCsparse = FALSE)
  }
  
  return(out)
}

#' Remove background contamination from count matrix
#' (refactored for maintainability, API unchanged)
adjustCounts = function(sc,clusters=NULL,method=c('subtraction','soupOnly','multinomial'),roundToInt=FALSE,verbose=1,tol=1e-3,pCut=0.01,...){
  method = match.arg(method)
  validate_soup_channel(sc, require_soup_profile = TRUE, require_rho = TRUE)
  if(!is.logical(roundToInt)) stop("roundToInt must be logical (TRUE/FALSE)")
  if(!is.numeric(verbose) || verbose < 0) stop("verbose must be a non-negative integer")
  if(!is.numeric(tol) || tol <= 0) stop("tol must be a positive number")
  if(!is.numeric(pCut) || pCut <= 0 || pCut >= 1) stop("pCut must be between 0 and 1")
  if(is.null(clusters)){
    if('clusters' %in% colnames(sc$metaData)){
      clusters = setNames(as.character(sc$metaData$clusters),rownames(sc$metaData))
    }else{
      warning("Clustering data not found.  Adjusting counts at cell level.  You will almost certainly get better results if you cluster data first.")
      clusters=FALSE
    }
  }
  if(clusters[1]!=FALSE){
    out = .adjust_counts_by_cluster(sc, clusters, method, verbose, tol, pCut, ...)
  }else{
    if(method=='multinomial'){
      if(verbose>1) message("Initialising with subtraction method.")
      fitInit = sc$toc - adjustCounts(sc,clusters=FALSE,method='subtraction',roundToInt=TRUE)
      ps = sc$soupProfile$est
      out = .fit_multinomial_counts(sc, fitInit, ps, verbose)
    }else if(method=='soupOnly'){
      out = .adjust_counts_soup_only(sc, pCut, verbose)
    }else if(method=='subtraction'){
      out = .adjust_counts_subtraction(sc)
    }else{
      stop("Internal error: Unknown method '", method, "' after validation. ",
           "This should not happen. Please report this bug.")
    }
  }
  if(roundToInt){
    if(verbose>1) message("Rounding to integers.")
    out@x = floor(out@x)+rbinom(length(out@x),1,out@x-floor(out@x))
  }
  return(out)
}

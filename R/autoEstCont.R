#' Automatically calculate the contamination fraction
#'
#' The idea of this method is that genes that are highly expressed in the soup and are marker genes for some population can be used to estimate the background contamination.  Marker genes are identified using the tfidf method (see \code{\link{quickMarkers}}).  The contamination fraction is then calculated at the cluster level for each of these genes and clusters are then aggressively pruned to remove those that give implausible estimates.
#'
#' This set of marker genes is filtered to include only those with tf-idf value greater than \code{tfidfMin}.  A higher tf-idf value implies a more specific marker.  Specifically a cut-off t implies that a marker gene has the property that geneFreqGlobal < exp(-t/geneFreqInClust).  See \code{\link{quickMarkers}}.  It may be necessary to decrease this value for data sets with few good markers.
#'
#' This set of marker genes is filtered down to include only the genes that are highly expressed in the soup, controlled by the \code{soupQuantile} parameter.  Genes highly expressed in the soup provide a more precise estimate of the contamination fraction.
#'
#' The pruning of implausible clusters is based on a call to \code{\link{estimateNonExpressingCells}}.  The parameters \code{maximumContamination=max(contaminationRange)} and \code{rhoMaxFDR} are passed to this function.  The defaults set here are calibrated to aggressively prune anything that has even the weakest of evidence that it is genuinely expressed. 
#'
#' For each cluster/gene pair the posterior distribution of the contamination fraction is calculated (based on gamma prior, controlled by \code{priorRho} and \code{priorRhoStdDev}).  These posterior distributions are aggregated to produce a final estimate of the contamination fraction. The logic behind this is that estimates from clusters that truly estimate the contamination fraction will cluster around the true value, while erroneous estimates will be spread out across the range (0,1) without a 'preferred value'.  The most probable value of the contamination fraction is then taken as the final global contamination fraction.
#'
#' @note
#' This function assumes that the channel contains multiple distinct cell types with different marker genes.  If you try and run it on a channel with very homogenous cells (e.g. a cell line, flow-sorted cells), you will likely get a warning, an error, and/or an extremely high contamination estimate.  In such circumstances your best option is usually to manually set the contamination to something reasonable.
#' 
#' @export
#' @param sc The SoupChannel object.
#' @param topMarkers A data.frame giving marker genes.  Must be sorted by decreasing specificity of marker and include a column 'gene' that contains the gene name.  If set to NULL, markers are estimated using \code{\link{quickMarkers}}.
#' @param tfidfMin Minimum value of tfidf to accept for a marker gene.
#' @param soupQuantile Only use genes that are at or above this expression quantile in the soup.  This prevents inaccurate estimates due to using genes with poorly constrained contribution to the background.
#' @param maxMarkers If we have heaps of good markers, keep only the best \code{maxMarkers} of them.
#' @param contaminationRange Vector of length 2 that constrains the contamination fraction to lie within this range.  Must be between 0 and 1.  The high end of this range is passed to \code{\link{estimateNonExpressingCells}} as \code{maximumContamination}.
#' @param rhoMaxFDR False discovery rate passed to \code{\link{estimateNonExpressingCells}}, to test if rho is less than \code{maximumContamination}.
#' @param priorRho Mode of gamma distribution prior on contamination fraction.
#' @param priorRhoStdDev Standard deviation of gamma distribution prior on contamination fraction.
#' @param doPlot Create a plot showing the density of estimates?
#' @param forceAccept Passed to \code{\link{setContaminationFraction}}.  Should we allow very high contamination fractions to be used.
#' @param verbose Be verbose?
#' @seealso quickMarkers
#' @return A modified SoupChannel object where the global contamination rate has been set.  Information about the estimation is also stored in the slot \code{fit}
#' @examples
#' #Use less specific markers
#' scToy = autoEstCont(scToy,tfidfMin=0.8)
#' #Allow large contamination fractions to be allocated
#' scToy = autoEstCont(scToy,forceAccept=TRUE)
#' #Be quiet
#' scToy = autoEstCont(scToy,verbose=FALSE,doPlot=FALSE)
#' @importFrom stats dgamma qgamma 
#' @importFrom graphics abline lines legend plot

# Internal helper: select and filter marker genes
.select_marker_genes <- function(sc, topMarkers, tfidfMin, soupQuantile, maxMarkers, verbose) {
  soupProf = sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ]
  soupMin = quantile(soupProf$est, soupQuantile)
  if (is.null(topMarkers)) {
    # Check if clusters are available, if not create a single cluster
    if(!"clusters" %in% colnames(sc$metaData)) {
      if(verbose) {
        message("No clustering information found. Creating single cluster for all cells.")
        message("For better results, consider providing clustering information using setClusters().")
      }
      # Create a single cluster for all cells
      sc$metaData$clusters <- rep("Cluster1", nrow(sc$metaData))
    }
    mrks = quickMarkers(sc$toc, sc$metaData$clusters, N = Inf)
    mrks = mrks[order(mrks$gene, -mrks$tfidf), ]
    mrks = mrks[!duplicated(mrks$gene), ]
    mrks = mrks[order(-mrks$tfidf), ]
    mrks = mrks[mrks$tfidf > tfidfMin, ]
  } else {
    mrks = topMarkers
  }
  tgts = rownames(soupProf)[soupProf$est > soupMin]
  filtPass = mrks[mrks$gene %in% tgts, ]
  tgts = head(filtPass$gene, n = maxMarkers)
  if (verbose)
    message(sprintf("%d genes passed tf-idf cut-off and %d soup quantile filter.  Taking the top %d.", nrow(mrks), nrow(filtPass), length(tgts)))
  list(mrks = mrks, tgts = tgts, soupProf = soupProf)
}

# Internal helper: optimized cluster-level rho estimation and posterior aggregation
.estimate_cluster_rho <- function(sc, ssc, tgts, mrks, soupProf, contaminationRange, rhoMaxFDR, priorRho, priorRhoStdDev, doPlot, verbose) {
  tmp = as.list(tgts)
  names(tmp) = tgts
  ute = estimateNonExpressingCells(sc, tmp, maximumContamination = max(contaminationRange), FDR = rhoMaxFDR)
  
  # Fix the logical subscript issue by properly handling the estimateNonExpressingCells output
  # ute has cell names as rows, we need to aggregate by cluster
  if(is.null(ssc$metaData$clusters) || length(ssc$metaData$clusters) == 1) {
    # Single cluster case - aggregate all cells
    cluster_cells = rownames(sc$metaData)
    if(length(cluster_cells) > 0) {
      ute_agg = matrix(colMeans(ute[cluster_cells, , drop = FALSE]), nrow = 1)
      rownames(ute_agg) = rownames(ssc$metaData)
      colnames(ute_agg) = colnames(ute)
      ute = t(ute_agg)
    } else {
      # Handle empty case
      ute = matrix(0, nrow = 1, ncol = ncol(ute))
      rownames(ute) = rownames(ssc$metaData)
      colnames(ute) = colnames(ute)
    }
  } else {
    # Multiple clusters case - aggregate by cluster
    cluster_names = rownames(ssc$metaData)
    ute_agg = matrix(0, nrow = length(cluster_names), ncol = ncol(ute))
    rownames(ute_agg) = cluster_names
    colnames(ute_agg) = colnames(ute)
    
    for(i in seq_along(cluster_names)) {
      cluster_name = cluster_names[i]
      cluster_cells = rownames(sc$metaData)[sc$metaData$clusters == cluster_name]
      if(length(cluster_cells) > 0) {
        ute_agg[i, ] = colMeans(ute[cluster_cells, , drop = FALSE])
      }
    }
    ute = t(ute_agg)
  }
  
  # Optimized expected counts calculation
  expCnts = outer(ssc$soupProfile$est, ssc$metaData$nUMIs)
  rownames(expCnts) = rownames(ssc$soupProfile)
  colnames(expCnts) = rownames(ssc$metaData)
  expCnts = expCnts[tgts, , drop = FALSE]
  obsCnts = ssc$toc[tgts, , drop = FALSE]
  
  # Vectorized p-value calculation
  pp = ppois(obsCnts, expCnts * max(contaminationRange), lower.tail = TRUE)
  qq = p.adjust(pp, method = 'BH')
  qq = matrix(qq, nrow = nrow(pp), ncol = ncol(pp), dimnames = dimnames(pp))
  
  # Vectorized rho estimation
  rhos = obsCnts / expCnts
  rhoIdx = t(apply(rhos, 1, function(e) order(order(e))))
  
  # More efficient data frame construction
  dd = data.frame(
    gene = rep(rownames(ute), ncol(ute)),
    passNonExp = as.vector(ute),
    rhoEst = as.vector(rhos),
    rhoIdx = as.vector(rhoIdx),
    obsCnt = as.vector(obsCnts),
    expCnt = as.vector(expCnts),
    isExpressedFDR = as.vector(qq),
    stringsAsFactors = FALSE
  )
  
  # Vectorized index matching
  dd$geneIdx = match(dd$gene, mrks$gene)
  dd$tfidf = mrks$tfidf[dd$geneIdx]
  dd$soupIdx = match(dd$gene, rownames(soupProf))
  dd$soupExp = soupProf$est[dd$soupIdx]
  dd$useEst = dd$passNonExp
  
  if (sum(dd$useEst) < 10)
    warning("Fewer than 10 independent estimates, rho estimation is likely to be unstable.  Consider reducing tfidfMin or increasing SoupMin.")
  if (verbose)
    message(sprintf("Using %d independent estimates of rho.", sum(dd$useEst)))
  
  # Optimized confidence interval calculation
  p.L = function(x, alpha) { if (x == 0) { 0 } else { qgamma(alpha, x) } }
  p.U = function(x, alpha) { qgamma(1 - alpha, x + 1) }
  alpha = 0.95
  alpha = (1 - alpha) / 2
  
  # Vectorized confidence interval calculation
  dd$rhoHigh = sapply(seq(nrow(dd)), function(e) p.U(dd$obsCnt[e], alpha) / dd$expCnt[e])
  dd$rhoLow = sapply(seq(nrow(dd)), function(e) p.L(dd$obsCnt[e], alpha) / dd$expCnt[e])
  
  # Optimized posterior calculation
  rhoProbes = seq(0, 1, .001)
  v2 = (priorRhoStdDev / priorRho) ** 2
  k = 1 + v2 ** -2 / 2 * (1 + sqrt(1 + 4 * v2))
  theta = priorRho / (k - 1)
  
  # Vectorized posterior calculation
  tmp = dd[dd$useEst, ]
  if(nrow(tmp) > 0) {
    post = sapply(rhoProbes, function(e) {
      mean(dgamma(e, k + tmp$obsCnt, scale = theta / (1 + theta * tmp$expCnt)))
    })
  } else {
    post = rep(0, length(rhoProbes))
  }
  
  xx = dgamma(rhoProbes, k, scale = theta)
  w = which(rhoProbes >= contaminationRange[1] & rhoProbes <= contaminationRange[2])
  
  if(length(w) > 0) {
    rhoEst = (rhoProbes[w])[which.max(post[w])]
    rhoFWHM = range((rhoProbes[w])[which(post[w] >= (max(post[w]) / 2))])
  } else {
    rhoEst = priorRho
    rhoFWHM = c(priorRho - priorRhoStdDev, priorRho + priorRhoStdDev)
  }
  
  contEst = rhoEst
  if (verbose)
    message(sprintf("Estimated global rho of %.2f", rhoEst))
  
  if (doPlot) {
    plot(rhoProbes, post, 'l',
         xlim = c(0, 1),
         ylim = c(0, max(c(xx, post))),
         frame.plot = FALSE,
         xlab = 'Contamination Fraction',
         ylab = 'Probability Density')
    lines(rhoProbes, xx, lty = 2)
    abline(v = rhoProbes[which.max(post)], col = 'red')
    legend(x = 'topright',
           legend = c(sprintf('prior rho %g(+/-%g)', priorRho, priorRhoStdDev),
                      sprintf('post rho %g(%g,%g)', rhoEst, rhoFWHM[1], rhoFWHM[2]),
                      'rho max'),
           lty = c(2, 1, 1),
           col = c('black', 'black', 'red'),
           bty = 'n')
  }
  
  list(dd = dd, post = post, rhoEst = rhoEst, rhoFWHM = rhoFWHM, markersUsed = mrks)
}

#' Automatically calculate the contamination fraction
#' (refactored for maintainability, API unchanged)
autoEstCont = function(sc,topMarkers=NULL,tfidfMin=1.0,soupQuantile=0.90,maxMarkers=100,contaminationRange=c(0.01,0.8),rhoMaxFDR=0.2,priorRho=0.05,priorRhoStdDev=0.10,doPlot=TRUE,forceAccept=FALSE,verbose=TRUE){
  validate_soup_channel(sc, require_soup_profile = TRUE, require_clusters = FALSE)
  if(length(contaminationRange) != 2 || any(contaminationRange < 0) || any(contaminationRange > 1) || contaminationRange[1] >= contaminationRange[2]) {
    stop("contaminationRange must be a vector of length 2 with values between 0 and 1, where first value < second value. ",
         "Got: [", paste(contaminationRange, collapse=", "), "]")
  }
  
  # Robust cluster aggregation for single or multiple clusters
  s = split(rownames(sc$metaData), sc$metaData$clusters)
  
  # Handle single cluster case properly
  if(length(s) == 1) {
    cluster_name <- names(s)[1]
    rs = rowSums(sc$toc[, s[[1]], drop = FALSE])
    tmp = matrix(rs, ncol = 1, nrow = nrow(sc$toc),
                dimnames = list(rownames(sc$toc), cluster_name))
  } else {
    tmp = do.call(cbind, lapply(s, function(e) {
      mat = sc$toc[, e, drop = FALSE]
      rs = rowSums(mat)
      if (is.null(dim(rs))) {
        rs = matrix(rs, ncol = 1, nrow = nrow(sc$toc),
                    dimnames = list(rownames(sc$toc), paste(e, collapse = "_")))
      }
      rs
    }))
    if (is.null(colnames(tmp))) colnames(tmp) <- names(s)
  }
  ssc = sc
  ssc$toc = tmp
  ssc$metaData = data.frame(nUMIs = colSums(tmp), row.names = colnames(tmp))
  marker_info = .select_marker_genes(ssc, topMarkers, tfidfMin, soupQuantile, maxMarkers, verbose)
  mrks = marker_info$mrks
  tgts = marker_info$tgts
  soupProf = marker_info$soupProf
  if(length(tgts) == 0){
    stop("No suitable marker genes found for contamination estimation. ",
         "Try: (1) reducing tfidfMin (currently ", tfidfMin, ") to accept less specific markers, ",
         "or (2) reducing soupQuantile (currently ", soupQuantile, ") to include lower-expressed genes, ",
         "or (3) if your data is homogeneous (e.g., cell line), manually set contamination with setContaminationFraction().")
  }
  if(length(tgts) < 10){
    warning("Fewer than 10 marker genes found.  Is this channel low complexity (see help)?  If not, consider reducing tfidfMin or soupQuantile")
  }
  est = .estimate_cluster_rho(sc, ssc, tgts, mrks, soupProf, contaminationRange, rhoMaxFDR, priorRho, priorRhoStdDev, doPlot, verbose)
  dd = est$dd
  post = est$post
  rhoEst = est$rhoEst
  rhoFWHM = est$rhoFWHM
  markersUsed = est$markersUsed
  sc$fit = list(dd = dd,
                priorRho = priorRho,
                priorRhoStdDev = priorRhoStdDev,
                posterior = post,
                rhoEst = rhoEst,
                rhoFWHM = rhoFWHM,
                markersUsed = markersUsed)
  sc = setContaminationFraction(sc, rhoEst, forceAccept = forceAccept)
  return(sc)
}

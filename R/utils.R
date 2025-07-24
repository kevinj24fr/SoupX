#' Expands soup counts calculated at the cluster level to the cell level
#'
#' Given a clustering of cells and soup counts calculated for each of those clusters,determines a most likely allocation of soup counts at the cell level.
#'
#' @export
#' @param clustSoupCnts Matrix of genes (rows) by clusters (columns) where counts are number of soup counts for that gene/cluster combination.
#' @param cellObsCnts Matrix of genes (rows) by cells (columns) giving the observed counts
#' @param clusters Mapping from cells to clusters.
#' @param cellWeights Weighting to give to each cell when distributing counts.  This would usually be set to the number of expected soup counts for each cell.
#' @param verbose Integer giving level of verbosity.  0 <- silence,1 <- Basic information,2 <- Very chatty,3 <- Debug.
#' @return A matrix of genes (rows) by cells (columns) giving the number of soup counts estimated for each cell.  Non-integer values possible.
#' @examples
#' # This function is typically called internally by adjustCounts()
#' # when using clustering information to improve correction accuracy
expandClusters <- function(clustSoupCnts,cellObsCnts,clusters,cellWeights,verbose=1){
  ws <- cellWeights
  #Do one cluster at a time
  if(verbose>0)
    message(sprintf("Expanding counts from %d clusters to %d cells.",ncol(clustSoupCnts),ncol(cellObsCnts)))
  
  # Pre-allocate result matrix for better performance
  result_matrix <- Matrix(0,nrow=nrow(cellObsCnts),ncol=ncol(cellObsCnts),sparse=TRUE)
  rownames(result_matrix) <- rownames(cellObsCnts)
  colnames(result_matrix) <- colnames(cellObsCnts)
  
  # Process clusters in parallel where possible
  cluster_names <- colnames(clustSoupCnts)
  
  for(j in seq_along(cluster_names)) {
    cluster_name <- cluster_names[j]
    if(verbose>1)
      message(sprintf("Expanding cluster %s",cluster_name))
    
    #Which cells
    wCells <- which(clusters==cluster_name)
    if(length(wCells) == 0) {
      if(verbose>1) message(sprintf("No cells found for cluster %s,skipping",cluster_name))
      next
    }
    
    #How should they be weighted
    ww <- ws[wCells]/sum(ws[wCells])
    
    #What is the limits
    lims <- cellObsCnts[,wCells,drop=FALSE]
    
    #And how many soup
    nSoup <- clustSoupCnts[,j]
    
    # Create sparse matrix directly without multiple conversions
    expCnts <- as(lims,"TsparseMatrix")
    
    #Most cases are easily dealt with.  In rough order of frequency.
    #1. No soup for gene - set to zero
    #2. All counts for gene are soup - set to lims
    #3. Not all counts are soup,but every entry is a 0 or 1 so no iteration needed.
    #4. Some iteration needed.
    
    #Deal with case 1 - vectorized operation
    zero_soup_genes <- which(nSoup==0)
    if(length(zero_soup_genes) > 0) {
      zero_indices <- which((expCnts@i+1) %in% zero_soup_genes)
      if(length(zero_indices) > 0) {
        expCnts@x[zero_indices] = 0
      }
    }
    
    #Case 2 is dealt with by construction
    
    #Get set of genes for cases 3 and 4
    wGenes <- which(nSoup>0 & nSoup<rowSums(lims))
    
    if(length(wGenes) > 0) {
      #And deal with them as appropriate. Save time by looking only at non-zero entries
      w <- which((expCnts@i+1) %in% wGenes)
      if(length(w) > 0) {
        w_split <- split(w,expCnts@i[w]+1)
        
        # Vectorized allocation for better performance
        tmp_results <- lapply(w_split,function(e) {
          gene_idx <- expCnts@i[e[1]]+1
          if(gene_idx <= length(nSoup)) {
            alloc(nSoup[gene_idx],expCnts@x[e],ww[expCnts@j[e]+1])
          } else {
            expCnts@x[e]  # Keep original values if gene index is out of bounds
          }
        })
        
        # Update values in one operation
        expCnts@x[unlist(w_split,use.names=FALSE)] = unlist(tmp_results,use.names=FALSE)
      }
    }
    
    # Assign to result matrix
    result_matrix[,wCells] <- expCnts
  }
  
  return(result_matrix)
}

#' Create Seurat style progress bar
#'
#' Creates progress bar that won't ruin log files and shows progress towards 100%.
#'
#' @param min Minimum value of parameter.
#' @param max Maximum value of parameter.
#' @param ... Passed to \code{\link{txtProgressBar}}
#' @return A txtProgressBar object to use updating progress.
initProgBar <- function(min,max,...){
  message('0%   10   20   30   40   50   60   70   80   90   100%')
  message('|----|----|----|----|----|----|----|----|----|----|')
  pb=txtProgressBar(min=min,max=max,style=1,width=51,char='*',...)
  return(pb)
}

#' Allocate values to "buckets" subject to weights and constraints
#'
#' Allocates \code{tgt} of something to \code{length(bucketLims)} different "buckets" subject to the constraint that each bucket has a maximum value of \code{bucketLims} that cannot be exceeded.  By default counts are distributed equally between buckets,but weights can be provided using \code{ws} to have the redistribution prefer certain buckets over others.
#'
#' @param tgt Value to distribute between buckets.
#' @param bucketLims The maximum value that each bucket can take.  Must be a vector of positive values.
#' @param ws Weights to be used for each bucket.  Default value makes all buckets equally likely.
#' @return A vector of the same length as \code{bucketLims} containing values distributed into buckets.
alloc <- function(tgt,bucketLims,ws=rep(1/length(bucketLims),length(bucketLims))){
  #Normalise weights
  ws <- ws/sum(ws)
  #Save time in line
  if(all(tgt*ws<=bucketLims))
    return(tgt*ws)
  #Need to order things in the order they'll be removed as the tgt increases
  o <- order(bucketLims/ws)
  w <- ws[o]
  y <- bucketLims[o]
  #The formula for number removed at entry i is
  #k_i <- \frac{y_i}{w_i} (1- \sum_j=0^{i-1} w_j) + \sum_j=0^{i-1} y_j
  cw <- cumsum(c(0,w[-length(w)]))
  cy <- cumsum(c(0,y[-length(y)]))
  k <- y/w* (1 - cw) + cy
  #Handle zero-weights appropriately
  k[w==0] = Inf
  #Everything that has k<=tgt will be set to y
  b <- (k<=tgt)
  #We then need to work out how many counts to distribute we have left over and distribute them according to re-normalised weights
  resid <- tgt-sum(y[b])
  w <- w/(1-sum(w[b]))
  out <- ifelse(b,y,resid*w)
  #Need to reverse sort
  return(out[order(o)])
}

#' Benchmark SoupX performance
#'
#' Measures execution time and memory usage for key SoupX operations to help users
#' understand performance characteristics and identify bottlenecks.
#'
#' @param sc A SoupChannel object to benchmark
#' @param operations Character vector of operations to benchmark. Options: 
#'   "quickMarkers","autoEstCont","adjustCounts","expandClusters","all"
#' @param iterations Number of iterations to run for each operation (default: 3)
#' @param verbose Whether to print detailed results
#' @return A data.frame with performance metrics for each operation
#' @examples
#' # Benchmark all operations
#' benchmark_soupx(scToy)
#' 
#' # Benchmark specific operations
#' benchmark_soupx(scToy,operations=c("quickMarkers","adjustCounts"))
#' @export
benchmark_soupx <- function(sc,operations="all",iterations=3,verbose=TRUE) {
  if(!inherits(sc,"SoupChannel")) {
    stop("sc must be a SoupChannel object",call. = FALSE)
  }
  
  if("all" %in% operations) {
    operations <- c("quickMarkers","autoEstCont","adjustCounts","expandClusters")
  }
  
  results <- data.frame(
    operation <- character(),
    mean_time <- numeric(),
    sd_time <- numeric(),
    min_time <- numeric(),
    max_time <- numeric(),
    memory_mb <- numeric(),
    stringsAsFactors <- FALSE
  )
  
  if(verbose) {
    message("Benchmarking SoupX operations...")
    message("Dataset size: ",nrow(sc$toc)," genes x ",ncol(sc$toc)," cells")
  }
  
  # Benchmark quickMarkers
  if("quickMarkers" %in% operations && "clusters" %in% colnames(sc$metaData)) {
    if(verbose) message("Benchmarking quickMarkers...")
    times <- numeric(iterations)
    memory_usage <- numeric(iterations)
    
    for(i in seq_len(iterations)) {
      gc() # Force garbage collection
      start_mem <- sum(gc()[,2])
      start_time <- Sys.time()
      
      mrks <- quickMarkers(sc$toc,sc$metaData$clusters,N=10)
      
      end_time <- Sys.time()
      end_mem <- sum(gc()[,2])
      
      times[i] <- as.numeric(difftime(end_time,start_time,units="secs"))
      memory_usage[i] <- (end_mem - start_mem) / 1024^2 # Convert to MB
    }
    
    results <- rbind(results,data.frame(
      operation <- "quickMarkers",
      mean_time <- mean(times),
      sd_time <- sd(times),
      min_time <- min(times),
      max_time <- max(times),
      memory_mb <- mean(memory_usage),
      stringsAsFactors <- FALSE
    ))
  }
  
  # Benchmark autoEstCont (requires clusters and soup profile)
  if("autoEstCont" %in% operations && 
     "clusters" %in% colnames(sc$metaData) && 
     !is.null(sc$soupProfile)) {
    if(verbose) message("Benchmarking autoEstCont...")
    times <- numeric(iterations)
    memory_usage <- numeric(iterations)
    
    for(i in seq_len(iterations)) {
      gc()
      start_mem <- sum(gc()[,2])
      start_time <- Sys.time()
      
      # Use tryCatch to handle potential errors gracefully
      result <- tryCatch({
        autoEstCont(sc,verbose=FALSE,doPlot=FALSE)
        TRUE
      },error <- function(e) {
        if(verbose) message("autoEstCont failed: ",e$message)
        FALSE
      })
      
      end_time <- Sys.time()
      end_mem <- sum(gc()[,2])
      
      if(result) {
        times[i] <- as.numeric(difftime(end_time,start_time,units="secs"))
        memory_usage[i] <- (end_mem - start_mem) / 1024^2
      } else {
        times[i] <- NA
        memory_usage[i] <- NA
      }
    }
    
    # Only add results if we have valid measurements
    valid_times <- times[!is.na(times)]
    if(length(valid_times) > 0) {
      results <- rbind(results,data.frame(
        operation <- "autoEstCont",
        mean_time <- mean(valid_times),
        sd_time <- sd(valid_times),
        min_time <- min(valid_times),
        max_time <- max(valid_times),
        memory_mb <- mean(memory_usage[!is.na(memory_usage)]),
        stringsAsFactors <- FALSE
      ))
    }
  }
  
  # Benchmark adjustCounts
  if("adjustCounts" %in% operations && 
     !is.null(sc$soupProfile) && 
     "rho" %in% colnames(sc$metaData)) {
    if(verbose) message("Benchmarking adjustCounts...")
    times <- numeric(iterations)
    memory_usage <- numeric(iterations)
    
    for(i in seq_len(iterations)) {
      gc()
      start_mem <- sum(gc()[,2])
      start_time <- Sys.time()
      
      result <- tryCatch({
        adjustCounts(sc,verbose=0)
        TRUE
      },error <- function(e) {
        if(verbose) message("adjustCounts failed: ",e$message)
        FALSE
      })
      
      end_time <- Sys.time()
      end_mem <- sum(gc()[,2])
      
      if(result) {
        times[i] <- as.numeric(difftime(end_time,start_time,units="secs"))
        memory_usage[i] <- (end_mem - start_mem) / 1024^2
      } else {
        times[i] <- NA
        memory_usage[i] <- NA
      }
    }
    
    valid_times <- times[!is.na(times)]
    if(length(valid_times) > 0) {
      results <- rbind(results,data.frame(
        operation <- "adjustCounts",
        mean_time <- mean(valid_times),
        sd_time <- sd(valid_times),
        min_time <- min(valid_times),
        max_time <- max(valid_times),
        memory_mb <- mean(memory_usage[!is.na(memory_usage)]),
        stringsAsFactors <- FALSE
      ))
    }
  }
  
  # Benchmark expandClusters (requires clusters)
  if("expandClusters" %in% operations && "clusters" %in% colnames(sc$metaData)) {
    if(verbose) message("Benchmarking expandClusters...")
    
    # Create dummy cluster-level data for testing
    clusters <- sc$metaData$clusters
    cluster_names <- unique(clusters)
    clustSoupCnts <- Matrix(0,nrow=nrow(sc$toc),ncol=length(cluster_names),sparse=TRUE)
    rownames(clustSoupCnts) <- rownames(sc$toc)
    colnames(clustSoupCnts) <- cluster_names
    
    # Fill with some dummy data
    for(i in seq_along(cluster_names)) {
      clustSoupCnts[,i] <- rowSums(sc$toc[,clusters == cluster_names[i],drop=FALSE]) * 0.1
    }
    
    times <- numeric(iterations)
    memory_usage <- numeric(iterations)
    
    for(i in seq_len(iterations)) {
      gc()
      start_mem <- sum(gc()[,2])
      start_time <- Sys.time()
      
      result <- tryCatch({
        expandClusters(clustSoupCnts,sc$toc,clusters,
                      sc$metaData$nUMIs * ifelse("rho" %in% colnames(sc$metaData),
                                                sc$metaData$rho,0.1),verbose=0)
        TRUE
      },error <- function(e) {
        if(verbose) message("expandClusters failed: ",e$message)
        FALSE
      })
      
      end_time <- Sys.time()
      end_mem <- sum(gc()[,2])
      
      if(result) {
        times[i] <- as.numeric(difftime(end_time,start_time,units="secs"))
        memory_usage[i] <- (end_mem - start_mem) / 1024^2
      } else {
        times[i] <- NA
        memory_usage[i] <- NA
      }
    }
    
    valid_times <- times[!is.na(times)]
    if(length(valid_times) > 0) {
      results <- rbind(results,data.frame(
        operation <- "expandClusters",
        mean_time <- mean(valid_times),
        sd_time <- sd(valid_times),
        min_time <- min(valid_times),
        max_time <- max(valid_times),
        memory_mb <- mean(memory_usage[!is.na(memory_usage)]),
        stringsAsFactors <- FALSE
      ))
    }
  }
  
  if(verbose) {
    message("\nPerformance Summary:")
    message("====================")
    for(i in seq_len(nrow(results))) {
      op <- results$operation[i]
      mean_t <- results$mean_time[i]
      mem <- results$memory_mb[i]
      message(sprintf("%-15s: %.3f ± %.3f sec,%.1f MB",op,mean_t,results$sd_time[i],mem))
    }
  }
  
  return(results)
}

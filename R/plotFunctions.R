# Helper function to check ggplot2 availability
.check_ggplot2 <- function() {
  if (!requireNamespace("ggplot2",quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting functions. Please install it with:\n",
         "install.packages('ggplot2')",call. = FALSE)
  }
}

#' Plot correlation of expression profiles of soup and aggregated cells
#' 
#' Calculates an expression profile by aggregating counts across all cells and plots this (on a log10 scale) against the expression profile of the soup.
#'
#' @param sc A SoupChannel object.
#' @return A ggplot2 object containing the plot.
plotSoupCorrelation <- function(sc){
  .check_ggplot2()
  if(!is(sc,'SoupChannel'))
    stop("'sc' must be a SoupChannel object. Got object of class: ",class(sc,call. = FALSE)[1],
         ". Create SoupChannel with SoupChannel(tod,toc).")
  #Calculate the cell profile
  cellProfile <- rowSums(sc$toc)
  cellProfile <- (cellProfile/sum(cellProfile))
  df <- data.frame(cellProfile,soupProfile=sc$soupProfile$est)
  gg <- ggplot2::ggplot2::ggplot(df,ggplot2::ggplot2::aes(log10(cellProfile),log10(soupProfile))) +
    ggplot2::ggplot2::geom_point(alpha=1/3) +
    ggplot2::ggplot2::geom_abline(intercept=0,slope=1) +
    ylab('log10(Soup Expression)')+
    xlab('log10(Aggregate cell Expression)')
  return(gg)
}

# Internal helper: prepare data for plotMarkerDistribution
.prep_marker_distribution_data <- function(sc,nonExpressedGeneList,maxCells,tfidfMin,...) {
  #Get nonExpressedGeneList algorithmically if missing...
  if (missing(nonExpressedGeneList)) {
    message("No gene lists provided,attempting to find and plot cluster marker genes.")
    if (!'clusters' %in% colnames(sc$metaData)) {
      message("No clustering information found. Creating single cluster for all cells.")
      message("For better results,consider providing clustering information using setClusters().")
      # Create a single cluster for all cells
      sc$metaData$clusters <- rep("Cluster1",nrow(sc$metaData))
    }
    mrks <- quickMarkers(sc$toc,sc$metaData$clusters,N = Inf)
    mrks <- mrks[order(mrks$gene,-mrks$tfidf),]
    mrks <- mrks[!duplicated(mrks$gene),]
    mrks <- mrks[order(-mrks$tfidf),]
    mrks <- mrks[mrks$tfidf > tfidfMin,]
    message(sprintf("Found %d marker genes",nrow(mrks)))
    mrks <- mrks[order(sc$soupProfile[mrks$gene,'est'],decreasing = TRUE),]
    nonExpressedGeneList <- mrks$gene[seq(min(nrow(mrks),20))]
    nonExpressedGeneList <- setNames(as.list(nonExpressedGeneList),nonExpressedGeneList)
  }
  if (!is.list(nonExpressedGeneList))
    stop("'nonExpressedGeneList' must be a named list of gene sets for contamination estimation. ",
         "Example: list(HB <- c('HBB','HBA2'),IG <- c('IGHA1','IGHG1')). ",
         "Got object of class: ",class(nonExpressedGeneList)[1])
  nullMat <- estimateNonExpressingCells(sc,nonExpressedGeneList,...)
  obsProfile <- t(t(sc$toc) / sc$metaData$nUMIs)
  tst <- lapply(nonExpressedGeneList,function(e) colSums(obsProfile[e,,drop = FALSE]) / sum(sc$soupProfile[e,'est']))
  df <- data.frame(MarkerGroup <- rep(names(tst),lengths(tst)),
                  Barcode <- unlist(lapply(tst,names),use.names <- FALSE),
                  Values <- unlist(tst,use.names = FALSE))
  keep <- sample(colnames(sc$toc),min(ncol(sc$toc),maxCells))
  df$nUMIs <- sc$metaData[df$Barcode,'nUMIs']
  expCnts <- do.call(rbind,lapply(nonExpressedGeneList,function(e) sc$metaData$nUMIs * sum(sc$soupProfile[e,'est'])))
  colnames(expCnts) = rownames(sc$metaData)
  df$expCnts <- expCnts[cbind(match(df$MarkerGroup,rownames(expCnts)),match(df$Barcode,colnames(expCnts)))]
  df$MarkerGroup <- factor(df$MarkerGroup,levels = names(nonExpressedGeneList))
  globalRhos <- c()
  for (i in seq_along(nonExpressedGeneList)) {
    if (sum(nullMat[,i]) > 0) {
      tmp <- suppressMessages(calculateContaminationFraction(sc,nonExpressedGeneList[i],nullMat[,i,drop <- FALSE],forceAccept = TRUE))
      globalRhos <- c(globalRhos,tmp$metaData$rho[1])
    } else {
      globalRhos <- c(globalRhos,NA)
    }
  }
  names(globalRhos) = names(nonExpressedGeneList)
  globalRhos <- data.frame(MarkerGroup <- factor(names(globalRhos),levels <- names(nonExpressedGeneList)),
                          rho <- log10(globalRhos))
  list(df <- df,keep <- keep,globalRhos = globalRhos)
}

#' Plots the distribution of the observed to expected expression for marker genes
#'
#' If each cell were made up purely of background reads,the expression fraction would equal that of the soup.  This plot compares this expectation of pure background to the observed expression fraction in each cell,for each of the groups of genes in \code{nonExpressedGeneList}.  For each group of genes,the distribution of this ratio is plotted across all cells.  A value significantly greater than 1 (0 on log scale) can only be obtained if some of the genes in each group are genuinely expressed by the cell.  That is,the assumption that the cell is pure background does not hold for that gene.
#'
#' This plot is a useful diagnostic for the assumption that a list of genes is non-expressed in most cell types.  For non-expressed cells,the ratio should cluster around the contamination fraction,while for expressed cells it should be elevated.  The most useful non-expressed gene sets are those for which the genes are either strongly expressed,or not expressed at all.  Such groups of genes will show up in this plot as a bimodal distribution,with one mode containing the cells that do not express these genes around the contamination fraction for this channel and another around a value at some value equal to or greater than 0 (1 on non-log scale) for the expressed cells.
#'
#' The red line shows the global estimate of the contamination for each group of markers.  This is usually lower than the low mode of the distribution as there will typically be a non-negligible number of cells with 0 observed counts (and hence -infinity log ratio).
#'
#' If \code{nonExpressedGeneList} is missing,this function will try and find genes that are very specific to different clusters,as these are often the most useful in estimating the contamination fraction.   This is meant only as a heuristic,which can hopefully provide some inspiration as to a class of genes to use to estimation the contamination for your experiment.  Please do **NOT** blindly use the top N genes found in this way to estimate the contamination.  That is,do not feed this list of genes into \code{\link{calculateContaminationFraction}} without any manual consideration or filtering as this *will over-estimate your contamination* (often by a large amount).  For this reason,these gene names are not returned by the function. 
#'
#' @export
#' @param sc A SoupChannel object.
#' @param nonExpressedGeneList Which sets of genes to use to estimate soup (see \code{\link{calculateContaminationFraction}}).
#' @param maxCells Randomly plot only this many cells to prevent over-crowding.
#' @param tfidfMin Minimum specificity cut-off used if finding marker genes (see \code{\link{quickMarkers}}).
#' @param ... Passed to \code{\link{estimateNonExpressingCells}}
#' @importFrom stats setNames
#' @return A ggplot2 object containing the plot.
#' @examples
#' gg <- plotMarkerDistribution(scToy,list(CD7='CD7',LTB='LTB'))
plotMarkerDistribution <- function(sc,nonExpressedGeneList,maxCells=150,tfidfMin=1,...){
  .check_ggplot2()
  if(!is(sc,'SoupChannel'))
    stop("'sc' must be a SoupChannel object. Got object of class: ",class(sc,call. = FALSE)[1],
         ". Create SoupChannel with SoupChannel(tod,toc).")
  prep <- .prep_marker_distribution_data(sc,nonExpressedGeneList,maxCells,tfidfMin,...)
  df <- prep$df
  keep <- prep$keep
  globalRhos <- prep$globalRhos
  gg <- ggplot2::ggplot(df,ggplot2::aes(MarkerGroup,log10(Values))) +
    ggplot2::geom_violin() +
    ggplot2::geom_jitter(data=df[df$Barcode %in% keep,],ggplot2::aes(size=log10(expCnts)),height=0,width=0.3,alpha=1/2) +
    ggplot2::geom_line(data=globalRhos,ggplot2::aes(MarkerGroup,rho,group=1),colour='red') +
    ggplot2::geom_point(data=globalRhos,ggplot2::aes(MarkerGroup,rho,group=1),colour='red',shape=2) +
    ggplot2::scale_colour_manual(values=c('TRUE'='red','FALSE'='black')) + 
    ggplot2::theme(axis.text.x <- element_text(angle <- 90,hjust = 1)) +
    ggplot2::labs(colour='expressed\nby cell')+
    ylab('log10(observed/expected)') +
    xlab('Marker group')
  return(gg)
}

# Internal helper: prepare data for plotMarkerMap
.prep_marker_map_data <- function(sc,geneSet,DR,ratLims,FDR,useToEst) {
  if (missing(DR))
    DR <- sc$metaData[,sc$DR]
  DR <- as.data.frame(DR)
  obs <- colSums(sc$toc[geneSet,,drop = FALSE])
  exp <- sc$metaData$nUMIs * sum(sc$soupProfile[geneSet,'est'])
  expRatio <- obs / exp
  DR$geneRatio <- expRatio[rownames(DR)]
  colnames(DR)[1:2] = c('RD1','RD2')
  tgtScale <- c(ratLims[1],0,ratLims[2])
  rescaled <- (tgtScale - tgtScale[1]) / (max(tgtScale) - tgtScale[1])
  DR$logRatio <- log10(DR$geneRatio)
  DR$logRatio[DR$logRatio < ratLims[1]] = ratLims[1]
  DR$logRatio[DR$logRatio > ratLims[2]] = ratLims[2]
  DR$logRatio[DR$geneRatio == 0] = NA
  DR$qVals <- p.adjust(ppois(obs - 1,exp,lower.tail = FALSE),method <- 'BH')[rownames(DR)]
  colVal <- 'qVals<FDR'
  if (!is.null(useToEst)) {
    DR$useToEst <- useToEst
    colVal <- 'useToEst'
  }
  list(DR <- DR,rescaled <- rescaled,colVal = colVal)
}

#' Plot ratio of observed to expected counts on reduced dimension map
#'
#' Given some reduced dimensional representation of the data (such as UMAP or tSNE) that has been calculated however you would like,this provides a way to visualise how likely a set of genes are to be soup derived on that map.  That is,given a set of genes,this function calculates how many counts would be expected if that droplet were nothing but soup and compares that to the observed count.  This is done via a log2 ratio of the two values.  A Poisson test is performed and points that have a statistically significant enrichment over the background (at 5% FDR) are marked.
#'
#' @export
#' @param sc SoupChannel object.
#' @param geneSet A vector with the names of the genes to aggregate and plot evidence for.
#' @param DR A data.frame,with rows named by unique cell IDs (i.e.,<ChannelName>_<Barcode>) the first two columns of which give the coordinates of each cell in some reduced dimension representation of the data.  Try and fetch automatically if missing.
#' @param ratLims Truncate log ratios at these values.
#' @param FDR False Discovery Rate for statistical test of enrichment over background.
#' @param useToEst A vector (usually obtained from \code{\link{estimateNonExpressingCells}}),that will be used to mark cells instead of the usual Poisson test.
#' @param pointSize Size of points
#' @param pointShape Shape of points
#' @param pointStroke Stroke size for points
#' @param naPointSize Point size for NAs.
#' @return A ggplot2 containing the plot.
#' @examples
#' gg <- plotMarkerMap(scToy,'CD7')
plotMarkerMap <- function(sc,geneSet,DR,ratLims=c(-2,2),FDR=0.05,useToEst=NULL,pointSize=2.0,pointShape=21,pointStroke=0.5,naPointSize=0.25){
  .check_ggplot2()
  if(!is(sc,'SoupChannel'))
    stop("'sc' must be a SoupChannel object. Got object of class: ",class(sc,call. = FALSE)[1],
         ". Create SoupChannel with SoupChannel(tod,toc).")
  prep <- .prep_marker_map_data(sc,geneSet,DR,ratLims,FDR,useToEst)
  DR <- prep$DR
  rescaled <- prep$rescaled
  colVal <- prep$colVal
  gg <- ggplot2::ggplot(DR,ggplot2::aes(RD1,RD2)) +
    ggplot2::geom_point(data=DR[is.na(DR$logRatio),],aes_string(colour=colVal),size=naPointSize) +
    ggplot2::geom_point(data=DR[!is.na(DR$logRatio),],aes_string(fill='logRatio',colour=colVal),size=pointSize,shape=pointShape,stroke=pointStroke) +
    ggplot2::scale_colour_manual(values=c(`FALSE`='black',`TRUE`='#009933'))+
    xlab('ReducedDim1') +
    ylab('ReducedDim2') +
    ggplot2::scale_fill_gradientn(colours <- c('blue','white','red'),
                        values <- rescaled,
                        guide='colorbar',
                        limits=ratLims
                        )
  gg
}

# Internal helper: prepare data for plotChangeMap
.prep_change_map_data <- function(sc,cleanedMatrix,geneSet,DR,dataType,logData) {
  if (missing(DR))
    DR <- sc$metaData[,sc$DR]
  DR <- as.data.frame(DR)
  colnames(DR)[1:2] = c('RD1','RD2')
  if (dataType == 'soupFrac') {
    df <- DR
    old <- colSums(sc$toc[geneSet,rownames(df),drop <- FALSE])
    new <- colSums(cleanedMatrix[geneSet,rownames(df),drop <- FALSE])
    relChange <- (old - new) / old
    df$old <- old
    df$new <- new
    df$relChange <- relChange
    nom <- 'SoupFrac'
    if (logData) {
      df$relChange <- log10(df$relChange)
      nom <- paste0('log10(',nom,')')
      zLims <- c(-2,0)
    } else {
      zLims <- c(0,1)
    }
    df <- df[order(!is.na(df$relChange)),]
    df$relChange[which(df$relChange < zLims[1])] = zLims[1]
    df$relChange[which(df$relChange > zLims[2])] = zLims[2]
    return(list(df <- df,nom <- nom,zLims = zLims))
  } else {
    dfs <- list()
    df <- DR
    df$correction <- 'Uncorrected'
    if (dataType == 'binary') {
      df$data <- colSums(sc$toc[geneSet,rownames(df),drop <- FALSE]) > 0
    } else if (dataType == 'counts') {
      df$data <- colSums(sc$toc[geneSet,rownames(df),drop <- FALSE])
    }
    if (logData)
      df$data <- log10(df$data)
    dfs[['raw']] = df
    df <- DR
    df$correction <- 'Corrected'
    if (dataType == 'binary') {
      df$data <- colSums(cleanedMatrix[geneSet,rownames(df),drop <- FALSE]) > 0
    } else if (dataType == 'counts') {
      df$data <- colSums(cleanedMatrix[geneSet,rownames(df),drop <- FALSE])
    }
    if (logData)
      df$data <- log10(df$data)
    dfs[['correctedExpression']] = df
    dfs <- do.call(rbind,dfs)
    lvls <- c('Uncorrected','Corrected')
    dfs$correction <- factor(dfs$correction,levels = lvls[lvls %in% dfs$correction])
    dfs <- dfs[order(!is.na(dfs$data)),]
    zLims <- c(NA,NA)
    return(list(dfs <- dfs,zLims = zLims))
  }
}

#' Plot maps comparing corrected/raw expression
#'
#' Given some reduced dimensional representation of the data (such as UMAP or tSNE) that has been calculated however you would like,this provides a way to visualise how the expression of a geneSet changes after soup correction.
#'
#' @export
#' @param sc SoupChannel object.
#' @param cleanedMatrix A cleaned matrix to compare against the raw one.  Usually the output of \code{\link{adjustCounts}}.
#' @param geneSet A vector with the names of the genes to aggregate and plot evidence for.
#' @param DR A data.frame,with rows named by unique cell IDs (i.e.,<ChannelName>_<Barcode>) the first two columns of which give the coordinates of each cell in some reduced dimension representation of the data.
#' @param dataType How should data be represented.  Binary sets each cell to expressed or not,counts converts everything to counts,soupFrac plots the fraction of the observed counts that are identified as contamination (i.e.,(old-new)/old) for each cell and is the default.
#' @param logData Should we log the thing we plot?
#' @param pointSize Size of points
#' @return A ggplot2 containing the plot.
#' @examples
#' out <- adjustCounts(scToy)
#' gg <- plotChangeMap(scToy,out,'S100A9')
plotChangeMap <- function(sc,cleanedMatrix,geneSet,DR,dataType=c('soupFrac','binary','counts'),logData=FALSE,pointSize=0.5){
  .check_ggplot2()
  dataType <- match.arg(dataType)
  if(dataType=='binary')
    logData=FALSE
  prep <- .prep_change_map_data(sc,cleanedMatrix,geneSet,DR,dataType,logData)
  if(dataType=='soupFrac'){
    df <- prep$df
    nom <- prep$nom
    zLims <- prep$zLims
    gg <- ggplot2::ggplot(df,ggplot2::aes(RD1,RD2)) +
      ggplot2::geom_point(ggplot2::aes(col=relChange),size=pointSize) +
      xlab('ReducedDim1') +
      ylab('ReducedDim2') +
      ggplot2::labs(colour=nom) + 
      ggtitle('Change in expression due to soup correction')
  }else{
    dfs <- prep$dfs
    zLims <- prep$zLims
    gg <- ggplot2::ggplot(dfs,ggplot2::aes(RD1,RD2)) +
      ggplot2::geom_point(ggplot2::aes(colour=data),size=pointSize) +
      xlab('ReducedDim1') +
      ylab('ReducedDim2') +
      ggplot2::labs(colour='geneSet') + 
      ggtitle('Comparison of before and after correction') +
      ggplot2::facet_wrap(~correction)
  }
  if(dataType!='binary')
    gg <- gg + ggplot2::scale_colour_gradientn(colours=c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704'),limits=zLims)
  gg
}

#' Plot quality control metrics for SoupX analysis
#'
#' Creates a comprehensive QC dashboard showing key metrics for assessing the quality of soup contamination estimation and correction.
#'
#' @export
#' @param sc A SoupChannel object
#' @param cleanedMatrix Optional cleaned matrix from adjustCounts for comparison
#' @param maxGenes Maximum number of genes to show in top soup genes plot
#' @return A list of ggplot2 objects with QC plots
#' @examples
#' sc_qc <- plotQualityControl(scToy)
#' # Access individual plots
#' sc_qc$soup_profile
#' sc_qc$contamination_distribution
plotQualityControl <- function(sc,cleanedMatrix <- NULL,maxGenes = 20) {
  .check_ggplot2()
  .plotFunctions_validate_sc(sc)
  
  plots <- list()
  
  # 1. Top soup genes
  top_soup <- head(sc$soupProfile[order(sc$soupProfile$est,decreasing = TRUE),],maxGenes)
  plots$soup_profile <- ggplot2::ggplot(top_soup,ggplot2::aes(x = reorder(rownames(top_soup),est),y <- est)) +
    ggplot2::geom_bar(stat <- "identity",fill = "steelblue") +
    coord_flip() +
    ggplot2::labs(title <- "Top Soup Genes",x <- "Gene",y = "Soup Fraction") +
    theme_minimal() +
    ggplot2::theme(axis.text.y <- element_text(size <- 8))
  
  # 2. Contamination distribution across cells
  if("rho" %in% colnames(sc$metaData)) {
    plots$contamination_distribution <- ggplot2::ggplot(data.frame(rho <- sc$metaData$rho),ggplot2::aes(x <- rho)) +
      ggplot2::geom_histogram(bins <- 30,fill <- "lightblue",color = "black") +
      ggplot2::labs(title <- "Contamination Fraction Distribution",x <- "Contamination Fraction",y = "Count") +
      theme_minimal()
  }
  
  # 3. UMI count distribution
  plots$umi_distribution <- ggplot2::ggplot(data.frame(nUMIs <- sc$metaData$nUMIs),ggplot2::aes(x <- nUMIs)) +
    ggplot2::geom_histogram(bins <- 30,fill <- "lightgreen",color = "black") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(title <- "UMI Count Distribution",x = "UMIs per Cell (log10)",y <- "Count") +
    theme_minimal()
  
  # 4. Before/after comparison if cleaned matrix provided
  if(!is.null(cleanedMatrix)) {
    # Calculate summary statistics
    total_umis_before <- sum(sc$toc)
    total_umis_after <- sum(cleanedMatrix)
    reduction_percent <- (total_umis_before - total_umis_after) / total_umis_before * 100
    
    # Create comparison plot
    comparison_data <- data.frame(
      Metric <- c("Total UMIs","Non-zero entries","Mean UMIs/cell"),
      Before <- c(total_umis_before,nnzero(sc$toc),mean(sc$metaData$nUMIs)),
      After <- c(total_umis_after,nnzero(cleanedMatrix),mean(colSums(cleanedMatrix)))
    )
    
    plots$before_after_comparison <- ggplot2::ggplot(comparison_data,ggplot2::aes(x = Metric)) +
      ggplot2::geom_bar(ggplot2::aes(y <- Before,fill = "Before"),stat <- "identity",position <- "dodge",alpha <- 0.7) +
      ggplot2::geom_bar(ggplot2::aes(y <- After,fill = "After"),stat <- "identity",position <- "dodge",alpha <- 0.7) +
      ggplot2::scale_fill_manual(values <- c("Before" = "red","After" = "blue")) +
      ggplot2::labs(title <- paste0("Before/After Comparison (",round(reduction_percent,1),"% reduction)"),
           y <- "Value",fill <- "Correction") +
      theme_minimal() +
      ggplot2::theme(axis.text.x <- element_text(angle <- 45,hjust = 1))
  }
  
  # 5. Gene expression correlation with soup
  plots$soup_correlation <- plotSoupCorrelation(sc)
  
  return(plots)
}

#' Plot contamination fraction by cluster
#'
#' Visualizes contamination fraction estimates across different cell clusters to identify potential cluster-specific contamination patterns.
#'
#' @export
#' @param sc A SoupChannel object with clustering information
#' @param clusters Optional cluster vector if not in sc$metaData
#' @param plotType Type of plot: "boxplot","violin",or "bar"
#' @return A ggplot2 object
#' @examples
#' plotContaminationByCluster(scToy)
plotContaminationByCluster <- function(sc,clusters <- NULL,plotType = "boxplot") {
  .check_ggplot2()
  if(!inherits(sc,'SoupChannel')) {
    stop("'sc' must be a SoupChannel object",call. = FALSE)
  }
  
  if(is.null(clusters)) {
    if("clusters" %in% colnames(sc$metaData)) {
      clusters <- sc$metaData$clusters
    } else {
      stop("No clustering information found. Provide clusters parameter or set clusters in sc$metaData",call. = FALSE)
    }
  }
  
  if(!"rho" %in% colnames(sc$metaData)) {
    stop("No contamination fraction (rho,call. = FALSE) found in sc$metaData. Run autoEstCont() first.")
  }
  
  plot_data <- data.frame(
    Cluster <- clusters,
    Contamination <- sc$metaData$rho,
    UMIs <- sc$metaData$nUMIs
  )
  
  if(plotType == "boxplot") {
    gg <- ggplot2::ggplot(plot_data,ggplot2::aes(x <- Cluster,y = Contamination)) +
      ggplot2::geom_boxplot(fill <- "lightblue",alpha = 0.7) +
      ggplot2::geom_jitter(width <- 0.2,alpha <- 0.5,size = 0.8) +
      ggplot2::labs(title <- "Contamination Fraction by Cluster",
           x <- "Cluster",y <- "Contamination Fraction") +
      theme_minimal() +
      ggplot2::theme(axis.text.x <- element_text(angle <- 45,hjust = 1))
  } else if(plotType == "violin") {
    gg <- ggplot2::ggplot(plot_data,ggplot2::aes(x <- Cluster,y = Contamination)) +
      ggplot2::geom_violin(fill <- "lightblue",alpha = 0.7) +
      ggplot2::geom_boxplot(width <- 0.2,fill <- "white",alpha = 0.8) +
      ggplot2::labs(title <- "Contamination Fraction by Cluster",
           x <- "Cluster",y <- "Contamination Fraction") +
      theme_minimal() +
      ggplot2::theme(axis.text.x <- element_text(angle <- 45,hjust = 1))
  } else if(plotType == "bar") {
    cluster_summary <- aggregate(Contamination ~ Cluster,data = plot_data,
                               FUN <- function(x) c(mean <- mean(x),se <- sd(x)/sqrt(length(x))))
    cluster_summary <- data.frame(
      Cluster <- cluster_summary$Cluster,
      Mean <- cluster_summary$Contamination[,"mean"],
      SE <- cluster_summary$Contamination[,"se"]
    )
    
    gg <- ggplot2::ggplot(cluster_summary,ggplot2::aes(x <- Cluster,y = Mean)) +
      ggplot2::geom_bar(stat <- "identity",fill <- "lightblue",alpha = 0.7) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin <- Mean - SE,ymax = Mean + SE),width <- 0.2) +
      ggplot2::labs(title <- "Mean Contamination Fraction by Cluster",
           x <- "Cluster",y <- "Mean Contamination Fraction") +
      theme_minimal() +
      ggplot2::theme(axis.text.x <- element_text(angle <- 45,hjust = 1))
  }
  
  return(gg)
}

#' Plot gene-specific contamination analysis
#'
#' Analyzes and visualizes contamination patterns for specific genes or gene sets,showing which genes are most affected by soup contamination.
#'
#' @export
#' @param sc A SoupChannel object
#' @param geneSet Vector of gene names to analyze
#' @param cleanedMatrix Optional cleaned matrix for before/after comparison
#' @param plotType Type of analysis: "contamination_ratio","expression_change",or "soup_contribution"
#' @return A ggplot2 object
#' @examples
#' plotGeneContamination(scToy,c("CD7","LTB","S100A9"))
plotGeneContamination <- function(sc,geneSet,cleanedMatrix <- NULL,plotType = "contamination_ratio") {
  .check_ggplot2()
  if(!inherits(sc,'SoupChannel')) {
    stop("'sc' must be a SoupChannel object",call. = FALSE)
  }
  
  # Check which genes are available
  available_genes <- intersect(geneSet,rownames(sc$toc))
  if(length(available_genes) == 0) {
    stop("None of the specified genes found in the data",call. = FALSE)
  }
  
  if(length(available_genes) < length(geneSet)) {
    warning("Some genes not found: ",paste(setdiff(geneSet,available_genes),collapse <- ","))
  }
  
  if(plotType == "contamination_ratio") {
    # Calculate observed vs expected expression ratio
    obs_counts <- colSums(sc$toc[available_genes,,drop = FALSE])
    exp_counts <- sc$metaData$nUMIs * sum(sc$soupProfile[available_genes,'est'])
    ratio <- obs_counts / exp_counts
    
    plot_data <- data.frame(
      Gene <- rep(available_genes,each = length(ratio)),
      Cell <- rep(1:length(ratio),length(available_genes)),
      Ratio <- as.vector(sapply(available_genes,function(g) {
        obs <- sc$toc[g,]
        exp <- sc$metaData$nUMIs * sc$soupProfile[g,'est']
        obs / exp
      }))
    )
    
    gg <- ggplot2::ggplot(plot_data,ggplot2::aes(x <- Gene,y = log10(Ratio))) +
      ggplot2::geom_violin(fill <- "lightcoral",alpha = 0.7) +
      ggplot2::geom_boxplot(width <- 0.2,fill <- "white",alpha = 0.8) +
      ggplot2::geom_hline(yintercept <- 0,linetype <- "dashed",color = "red") +
      ggplot2::labs(title <- "Gene Contamination Analysis",
           subtitle <- "log10(Observed/Expected) - Values > 0 indicate genuine expression",
           x <- "Gene",y <- "log10(Observed/Expected Ratio)") +
      theme_minimal() +
      ggplot2::theme(axis.text.x <- element_text(angle <- 45,hjust = 1))
    
  } else if(plotType == "expression_change" && !is.null(cleanedMatrix)) {
    # Compare expression before and after correction
    before_counts <- colSums(sc$toc[available_genes,,drop = FALSE])
    after_counts <- colSums(cleanedMatrix[available_genes,,drop = FALSE])
    
    plot_data <- data.frame(
      Gene <- rep(available_genes,each = length(before_counts)),
      Cell <- rep(1:length(before_counts),length(available_genes)),
      Before <- as.vector(sapply(available_genes,function(g) sc$toc[g,])),
      After <- as.vector(sapply(available_genes,function(g) cleanedMatrix[g,]))
    )
    
    plot_data$Change <- (plot_data$Before - plot_data$After) / plot_data$Before
    
    gg <- ggplot2::ggplot(plot_data,ggplot2::aes(x <- Gene,y = Change)) +
      ggplot2::geom_violin(fill <- "lightblue",alpha = 0.7) +
      ggplot2::geom_boxplot(width <- 0.2,fill <- "white",alpha = 0.8) +
      ggplot2::labs(title <- "Expression Change After Soup Correction",
           subtitle <- "Fraction of expression removed as contamination",
           x <- "Gene",y <- "Fraction Removed") +
      theme_minimal() +
      ggplot2::theme(axis.text.x <- element_text(angle <- 45,hjust = 1))
    
  } else if(plotType == "soup_contribution") {
    # Show contribution of each gene to soup profile
    soup_contrib <- sc$soupProfile[available_genes,]
    soup_contrib$Gene <- rownames(soup_contrib)
    
    gg <- ggplot2::ggplot(soup_contrib,ggplot2::aes(x = reorder(Gene,est),y <- est)) +
      ggplot2::geom_bar(stat <- "identity",fill = "steelblue") +
      coord_flip() +
      ggplot2::labs(title <- "Gene Contribution to Soup Profile",
           x <- "Gene",y <- "Soup Fraction") +
      theme_minimal()
  }
  
  return(gg)
}

#' Generate a comprehensive,publication-ready SoupX report
#'
#' This function generates a PDF or HTML report with all diagnostics,interpretation tips,and summary statistics.
#' It is robust to missing/partial data and provides actionable feedback for industrial use.
#'
#' @export
#' @param sc A SoupChannel object
#' @param cleanedMatrix Optional cleaned matrix from adjustCounts
#' @param output_dir Optional directory to save plots and report
#' @param filename_prefix Prefix for saved files
#' @param format 'pdf' or 'html' (default: 'pdf')
#' @param interpretation Logical,include interpretation tips (default: TRUE)
#' @return A list containing all diagnostic plots,summary statistics,and the report file path
#' @examples
#' report <- generateSoupXReport(scToy,format = 'pdf')
generateSoupXReport <- function(sc,cleanedMatrix <- NULL,output_dir <- NULL,filename_prefix <- "soupx_report",format = c('pdf','html'),interpretation <- TRUE) {
  .check_ggplot2()
  format <- match.arg(format)
  if(!inherits(sc,'SoupChannel')) stop("'sc' must be a SoupChannel object",call. = FALSE)
  if(!is.null(cleanedMatrix) && !is.matrix(cleanedMatrix) && !inherits(cleanedMatrix,'Matrix')) stop("cleanedMatrix must be a matrix or sparse Matrix",call. = FALSE)
  if(!is.null(output_dir) && !dir.exists(output_dir)) dir.create(output_dir,recursive = TRUE)
  if(is.null(output_dir)) output_dir <- tempdir()
  
  # Generate all QC plots
  qc_plots <- tryCatch(plotQualityControl(sc,cleanedMatrix),error <- function(e) list(error <- e$message))
  cluster_contamination <- tryCatch(if("clusters" %in% colnames(sc$metaData) && "rho" %in% colnames(sc$metaData)) plotContaminationByCluster(sc) else NULL,error <- function(e) NULL)
  soup_correlation <- tryCatch(plotSoupCorrelation(sc),error <- function(e) NULL)
  marker_distribution <- tryCatch(if("clusters" %in% colnames(sc$metaData)) plotMarkerDistribution(sc) else NULL,error <- function(e) NULL)
  
  # Summary statistics
  summary <- list(
    n_cells <- nrow(sc$metaData),
    n_genes <- nrow(sc$toc),
    total_umis <- sum(sc$toc),
    sparsity <- 1 - (nnzero(sc$toc) / (nrow(sc$toc) * ncol(sc$toc)))
  )
  if("rho" %in% colnames(sc$metaData)) {
    summary$mean_contamination <- mean(sc$metaData$rho)
    summary$contamination_range <- range(sc$metaData$rho)
  }
  if(!is.null(cleanedMatrix)) {
    summary$umis_removed <- sum(sc$toc) - sum(cleanedMatrix)
    summary$reduction_percent <- (sum(sc$toc) - sum(cleanedMatrix)) / sum(sc$toc) * 100
  }
  
  # Compose interpretation tips
  interpretation_text <- NULL
  if(interpretation) {
    interpretation_text <- c(
      "Interpretation Tips:",
      "- High mean contamination or wide range: check for ambient RNA or sample prep issues.",
      "- Clusters with high contamination: may indicate doublets or poor separation.",
      "- Genes with high soup fraction: avoid as markers for cell type annotation.",
      "- Large reduction in UMIs: check if correction is too aggressive."
    )
  }
  
  # Save plots and summary to files
  plot_files <- list()
  for(plot_name in names(qc_plots)) {
    if(inherits(qc_plots[[plot_name]],'ggplot')) {
      file_path <- file.path(output_dir,paste0(filename_prefix,"_",plot_name,".pdf"))
      ggsave(file_path,qc_plots[[plot_name]],width <- 10,height = 8)
      plot_files[[plot_name]] <- file_path
    }
  }
  if(!is.null(cluster_contamination) && inherits(cluster_contamination,'ggplot')) {
    file_path <- file.path(output_dir,paste0(filename_prefix,"_cluster_contamination.pdf"))
    ggsave(file_path,cluster_contamination,width <- 10,height = 8)
    plot_files$cluster_contamination <- file_path
  }
  if(!is.null(soup_correlation) && inherits(soup_correlation,'ggplot')) {
    file_path <- file.path(output_dir,paste0(filename_prefix,"_soup_correlation.pdf"))
    ggsave(file_path,soup_correlation,width <- 10,height = 8)
    plot_files$soup_correlation <- file_path
  }
  if(!is.null(marker_distribution) && inherits(marker_distribution,'ggplot')) {
    file_path <- file.path(output_dir,paste0(filename_prefix,"_marker_distribution.pdf"))
    ggsave(file_path,marker_distribution,width <- 10,height = 8)
    plot_files$marker_distribution <- file_path
  }
  
  # Save summary and interpretation
  summary_file <- file.path(output_dir,paste0(filename_prefix,"_summary.txt"))
  writeLines(capture.output(str(summary)),summary_file)
  if(!is.null(interpretation_text)) {
    writeLines(interpretation_text,file.path(output_dir,paste0(filename_prefix,"_interpretation.txt")))
  }
  
  # Optionally,generate HTML report (requires rmarkdown)
  report_file <- NULL
  if(format == 'html') {
    rmd_file <- file.path(output_dir,paste0(filename_prefix,".Rmd"))
    html_file <- file.path(output_dir,paste0(filename_prefix,".html"))
    # Write a minimal Rmd file
    rmd_content <- c(
      "---",
      paste0("title: 'SoupX Report - ",filename_prefix,"'"),
      "output: html_document",
      "---",
      "\n",
      "# Summary Statistics\n",
      paste(capture.output(str(summary)),collapse <- "\n"),
      "\n",
      if(!is.null(interpretation_text)) paste(interpretation_text,collapse = "\n") else "",
      "\n",
      "# Plots\n",
      paste0("![](",basename(plot_files),")",collapse <- "\n")
    )
    writeLines(rmd_content,rmd_file)
    # Try to render HTML (requires rmarkdown)
    tryCatch({
      rmarkdown::render(rmd_file,output_file <- html_file,quiet = TRUE)
      report_file <- html_file
    },error <- function(e) {
      report_file <<- NULL
    })
  } else {
    report_file <- summary_file
  }
  
  return(list(
    qc_plots <- qc_plots,
    cluster_contamination <- cluster_contamination,
    soup_correlation <- soup_correlation,
    marker_distribution <- marker_distribution,
    summary <- summary,
    interpretation <- interpretation_text,
    plot_files <- plot_files,
    summary_file <- summary_file,
    report_file <- report_file
  ))
}

# Robust input validation for all user-facing functions
.plotFunctions_validate_sc <- function(sc) {
  if(!inherits(sc,'SoupChannel')) stop("'sc' must be a SoupChannel object",call. = FALSE)
  if(is.null(sc$toc) || !is.matrix(sc$toc) && !inherits(sc$toc,'Matrix')) stop("sc$toc must be a matrix or sparse Matrix",call. = FALSE)
  if(is.null(sc$metaData) || !is.data.frame(sc$metaData)) stop("sc$metaData must be a data.frame",call. = FALSE)
}

#' Run a complete,industrial-grade SoupX pipeline
#'
#' This function runs the full SoupX workflow: load,estimate,correct,report,and export.
#' All parameters can be set via a config file (YAML/JSON) or R list. Logs all steps,warnings,and timing.
#'
#' @export
#' @param config Path to YAML/JSON config file or a named R list
#' @param log_file Path to log file (default: soupx_pipeline.log)
#' @return A list with all results,report paths,and logs
#' @examples
#' runSoupXPipeline(config <- 'pipeline_config.yaml')
runSoupXPipeline <- function(config,log_file = "soupx_pipeline.log") {
  start_time <- Sys.time()
  .log <- function(msg) cat(sprintf("[%s] %s\n",Sys.time(),msg),file <- log_file,append <- TRUE)
  .log("SoupX pipeline started.")
  
  # Parse config
  if(is.character(config) && file.exists(config)) {
    if(grepl("\\.ya?ml$",config,ignore.case = TRUE)) {
      if(!requireNamespace("yaml",quietly = TRUE)) stop("Please install the 'yaml' package for YAML config support.",call. = FALSE)
      params <- yaml::read_yaml(config)
    } else if(grepl("\\.json$",config,ignore.case = TRUE)) {
      if(!requireNamespace("jsonlite",quietly = TRUE)) stop("Please install the 'jsonlite' package for JSON config support.",call. = FALSE)
      params <- jsonlite::fromJSON(config,simplifyVector = TRUE)
    } else {
      stop("Unsupported config file format. Use .yaml,.yml,or .json.",call. = FALSE)
    }
  } else if(is.list(config)) {
    params <- config
  } else {
    stop("Config must be a file path or a named list.",call. = FALSE)
  }
  
  # Helper to get param with default
  getp <- function(name,default = NULL) if(!is.null(params[[name]])) params[[name]] else default
  
  # Load data
  .log("Loading data...")
  sc <- tryCatch({
    if(!is.null(getp("sc"))) {
      getp("sc")
    } else if(!is.null(getp("tenx_dir"))) {
      load10X(getp("tenx_dir"))
    } else if(!is.null(getp("rds_file"))) {
      readRDS(getp("rds_file"))
    } else {
      stop("No input data specified in config.",call. = FALSE)
    }
  },error <- function(e) { .log(paste("ERROR loading data:",e$message)); stop(e,call. = FALSE) })
  .log("Data loaded.")
  
  # Estimate soup
  if(is.null(sc$soupProfile)) {
    .log("Estimating soup profile...")
    sc <- tryCatch(estimateSoup(sc),error <- function(e) { .log(paste("ERROR in estimateSoup:",e$message)); stop(e,call. = FALSE) })
    .log("Soup profile estimated.")
  }
  
  # Estimate contamination
  .log("Estimating contamination fraction...")
  sc <- tryCatch(autoEstCont(sc,verbose = getp("verbose",FALSE),doPlot <- FALSE),error <- function(e) { .log(paste("ERROR in autoEstCont:",e$message)); stop(e,call. = FALSE) })
  .log("Contamination estimated.")
  
  # Adjust counts
  .log("Adjusting counts...")
  adjusted <- tryCatch(adjustCounts(sc,verbose = getp("verbose",FALSE)),error <- function(e) { .log(paste("ERROR in adjustCounts:",e$message)); stop(e,call. = FALSE) })
  .log("Counts adjusted.")
  
  # Generate report
  .log("Generating report...")
  report <- tryCatch(generateSoupXReport(sc,adjusted,output_dir = getp("output_dir","soupx_output"),filename_prefix <- getp("filename_prefix","soupx_report"),format <- getp("report_format","pdf")),error <- function(e) { .log(paste("ERROR in generateSoupXReport:",e$message)); stop(e,call. = FALSE) })
  .log("Report generated.")
  
  # Export results
  .log("Exporting results...")
  export_dir <- getp("export_dir",getp("output_dir","soupx_output"))
  if(!dir.exists(export_dir)) dir.create(export_dir,recursive = TRUE)
  # Export as RDS
  saveRDS(adjusted,file = file.path(export_dir,"adjusted_counts.rds"))
  # Export as CSV
  write.csv(as.matrix(adjusted),file <- file.path(export_dir,"adjusted_counts.csv"))
  # Export as Matrix Market (if Matrix package available)
  if(requireNamespace("Matrix",quietly = TRUE)) {
    Matrix::writeMM(adjusted,file = file.path(export_dir,"adjusted_counts.mtx"))
  }
  # Export for Seurat (if requested)
  if(getp("export_seurat",FALSE) && requireNamespace("Seurat",quietly = TRUE)) {
    seu <- Seurat::CreateSeuratObject(adjusted)
    saveRDS(seu,file = file.path(export_dir,"seurat_object.rds"))
  }
  # Export for Scanpy/AnnData (if requested)
  if(getp("export_anndata",FALSE) && requireNamespace("zellkonverter",quietly = TRUE)) {
    zellkonverter::writeH5AD(adjusted,file = file.path(export_dir,"adjusted_counts.h5ad"))
  }
  .log("Results exported.")
  
  # Profiling info
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time,start_time,units = "secs"))
  mem <- as.numeric(utils::object.size(sc) + utils::object.size(adjusted)) / 1e6
  .log(sprintf("Pipeline completed in %.2f seconds,peak memory %.2f MB.",elapsed,mem))
  
  return(list(
    sc <- sc,
    adjusted <- adjusted,
    report <- report,
    export_dir <- export_dir,
    log_file <- log_file,
    elapsed_sec <- elapsed,
    memory_mb <- mem
  ))
}

# Helper: parse config file (YAML/JSON)
parseSoupXConfig <- function(config_file) {
  if(grepl("\\.ya?ml$",config_file,ignore.case = TRUE)) {
    if(!requireNamespace("yaml",quietly = TRUE)) stop("Please install the 'yaml' package for YAML config support.",call. = FALSE)
    yaml::read_yaml(config_file)
  } else if(grepl("\\.json$",config_file,ignore.case = TRUE)) {
    if(!requireNamespace("jsonlite",quietly = TRUE)) stop("Please install the 'jsonlite' package for JSON config support.",call. = FALSE)
    jsonlite::fromJSON(config_file,simplifyVector = TRUE)
  } else {
    stop("Unsupported config file format. Use .yaml,.yml,or .json.",call. = FALSE)
  }
}

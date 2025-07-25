#' Load a collection of 10X data-sets
#'
#' Loads unfiltered 10X data from each data-set and identifies which droplets are cells using the cellranger defaults.
#'
#' @export
#' @param dataDir Top level cellranger output directory (the directory that contains the \code{raw_gene_bc_matrices} folder).
#' @param cellIDs Barcodes of droplets that contain cells.  If NULL,use the default cellranger set.
#' @param channelName The name of the channel to store.  If NULL set to either \code{names(dataDir)} or \code{dataDir} is no name is set.
#' @param readArgs A list of extra parameters passed to \code{Seurat::Read10X}.
#' @param includeFeatures If multiple feature types are present,keep only the types mentioned here and collapse to a single matrix.
#' @param verbose Be verbose?
#' @param ... Extra parameters passed to \code{SoupChannel} construction function.
#' @return A SoupChannel object containing the count tables for the 10X dataset.
#' @seealso SoupChannel estimateSoup
#' @examples
#' sc <- load10X(system.file('extdata','toyData',package='SoupX'))
#' @importFrom utils read.csv
load10X <- function(dataDir,cellIDs=NULL,channelName=NULL,readArgs=list(),includeFeatures=c('Gene Expression'),verbose=TRUE,...){
  # Check if Seurat is available
  if (!requireNamespace("Seurat",quietly = TRUE)) {
    stop("Package 'Seurat' is required for load10X function. Please install it with:\n",
         "install.packages('Seurat')",call. = FALSE)
  }
  #Work out which version we're dealing with
  isV3 <- dir.exists(file.path(dataDir,'raw_feature_bc_matrix'))
  isV7 <- dir.exists(file.path(dataDir,'analysis','clustering','gene_expression_graphclust'))
  isMulti <- dir.exists(file.path(dataDir,'analysis','clustering','gex'))
  tgt <- file.path(dataDir,
                  ifelse(isV3,'raw_feature_bc_matrix','raw_gene_bc_matrices'))
  #Add the reference genome for the non-V3 ones
  if(!isV3)
    tgt <- file.path(tgt,list.files(tgt))
  if(verbose)
    message(sprintf("Loading raw count data"))
  dat <- do.call(Seurat::Read10X,c(list(data.dir=tgt),readArgs))
  if(verbose)
    message(sprintf("Loading cell-only count data"))
  if(!is.null(cellIDs)){
    #This is now redundant as we require the version of Seurat that does not strip the suffix
    #Check we have the IDs
    if(!all(cellIDs %in% colnames(dat)))
      stop("Not all supplied cellIDs found in raw data. ",
           "Missing cellIDs: ",paste(head(setdiff(cellIDs,colnames(dat)),10),collapse=","),
           ". Check cellID format matches raw data column names.")
    datCells <- dat[,match(cellIDs,colnames(dat))]
  }else{
    #Work out which ones contain cells
    tgt <- file.path(dataDir,
                    ifelse(isV3,'filtered_feature_bc_matrix','filtered_gene_bc_matrices'))
    if(!isV3)
      tgt <- file.path(tgt,list.files(tgt))
    datCells <- do.call(Seurat::Read10X,c(list(data.dir=tgt),readArgs))
    #If it's a list of multiple types,have to decide what to include and collapse to one matrix.
    if(is.list(dat)){
      dat <- do.call(rbind,dat[includeFeatures])
      datCells <- do.call(rbind,datCells[includeFeatures])
    }
  }
  if(verbose)
    message(sprintf("Loading extra analysis data where available"))
  #Get the cluster annotation if available
  mDat <- NULL
  #What needs to be added to make V7 directory structure work
  v7Prefix=ifelse(isV7,'gene_expression_','')
  tgt <- ifelse(isMulti && !isV7,
               file.path(dataDir,'analysis','clustering','gex','graphclust','clusters.csv'),
               file.path(dataDir,'analysis','clustering',paste0(v7Prefix,'graphclust'),'clusters.csv')
               )
  if(file.exists(tgt)){
    clusters <- read.csv(tgt)
    mDat <- data.frame(clusters=clusters$Cluster,row.names=clusters$Barcode)
  }
  #Add fine grained clusters too if present
  tgt <- ifelse(isMulti && !isV7,
               file.path(dataDir,'analysis','clustering','gex','kmeans_10_clusters','clusters.csv'),
               file.path(dataDir,'analysis','clustering',paste0(v7Prefix,'kmeans_10_clusters'),'clusters.csv')
               )
  if(file.exists(tgt)){
    clusters <- read.csv(tgt)
    mDat$clustersFine <- clusters$Cluster
  }
  #Get tSNE if available and point to it
  tgt <- ifelse(isMulti && !isV7,
               file.path(dataDir,'analysis','dimensionality_reduction','gex','tsne_projection.csv'),
               file.path(dataDir,'analysis','tsne',paste0(v7Prefix,'2_components'),'projection.csv')
               )
  if(file.exists(tgt)){
    tsne <- read.csv(tgt)
    if(is.null(mDat)){
      mDat <- data.frame(tSNE1=tsne$TSNE.1,tSNE2=tsne$TSNE.2,row.names=tsne$Barcode)
    }else{
      mDat$tSNE1 <- tsne$TSNE.1[match(rownames(mDat),tsne$Barcode)]
      mDat$tSNE2 <- tsne$TSNE.2[match(rownames(mDat),tsne$Barcode)]
    }
    DR <- c('tSNE1','tSNE2')
  }else{
    DR=NULL
  }
  #Ensure rownames of metadata match column names of counts
  if(!is.null(mDat) && any(rownames(mDat)!=colnames(datCells))){
    rownames(mDat) = gsub('-1$','',rownames(mDat))
    if(any(rownames(mDat)!=colnames(datCells)))
      stop("Error matching metadata to cell names. ",
           "Metadata rownames must match cell names in count matrix. ",
           "Check that cell barcodes are consistent between data files.")
  }
  #Get a name for the channel
  if(is.null(channelName))
    channelName <- ifelse(is.null(names(dataDir)),dataDir,names(dataDir))
  channel <- SoupChannel(tod = dat,
                        toc = datCells,
                        metaData = mDat,
                        channelName = channelName,
                        dataDir = dataDir,
                        dataType='10X',
                        isV3=isV3,
                        DR=DR,
                        ...
                        )
  return(channel)
}

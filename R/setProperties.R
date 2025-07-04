#' Set soup profile
#'
#' Manually sets or updates the soup profile for a SoupChannel object. 
#'
#' @export
#' @param sc A SoupChannel object.
#' @param soupProfile A data.frame with columns \code{est} containing the fraction of soup for each gene, \code{counts} containing the total counts for each gene and with row names corresponding to the row names of \code{sc$toc}.
#' @return An updated SoupChannel object with the soup profile set.
#' @examples
#' #Suppose only table of counts is available
#' toc = Seurat::Read10X(system.file('extdata','toyData','filtered_gene_bc_matrices','GRCh38',
#'                                   package='SoupX'))
#' #Suppress calculating soup profile automatically
#' sc = SoupChannel(toc,toc,calcSoupProfile=FALSE)
#' #And add manually
#' rowSums = Matrix::rowSums
#' soupProf = data.frame(row.names = rownames(toc),est=rowSums(toc)/sum(toc),counts=rowSums(toc))
#' sc = setSoupProfile(sc,soupProf)
setSoupProfile = function(sc,soupProfile){
  # Validate inputs
  validate_soup_channel(sc)
  validate_soup_profile(soupProfile, rownames(sc$toc))
  
  if(!all(rownames(soupProfile) %in% rownames(sc$toc))){
    stop("soupProfile missing genes found in count matrix. ",
         "Missing genes: ", paste(head(setdiff(rownames(sc$toc), rownames(soupProfile)), 10), collapse=", "))
  }else{
    sc$soupProfile = soupProfile[rownames(sc$toc),]
  }
  return(sc)
}

#' Sets clustering for SoupChannel
#' 
#' Adds or updates clustering information to meta-data table in SoupChannel object.
#'
#' @export
#' @param sc A SoupChannel object.
#' @param clusters A named vector, where entries are the cluster IDs and names are cellIDs.  If no names are provided, the order is assumed to match the order in \code{sc$metaData}.
#' @return An updated SoupChannel object with clustering information stored.
#' @examples
#' sc = load10X(system.file('extdata','toyData',package='SoupX'))
#' mDat = read.table(system.file('extdata','toyData','metaData.tsv',package='SoupX'),sep='\t')
#' sc = setClusters(sc,mDat$res.1)
setClusters = function(sc,clusters){
  # Validate inputs
  validate_soup_channel(sc)
  validate_clusters(clusters, rownames(sc$metaData))
  
  if(!all(colnames(sc$toc) %in% names(clusters))){
    if(length(clusters)==nrow(sc$metaData)){
      #Ensure the thing we're setting is not a factor
      sc$metaData$clusters = as.character(clusters)
    }else{
      stop("Invalid cluster specification: 'clusters' must be either a named vector with all cell names, ",
           "or an unnamed vector with length equal to number of cells (", nrow(sc$metaData), "). ",
           "Got length: ", length(clusters))
    }
  }else{
    sc$metaData$clusters = as.character(clusters[rownames(sc$metaData)])
  }
  return(sc)
}

#' Manually set contamination fraction
#'
#' Manually specify the contamination fraction.
#'
#' @export
#' @param sc A SoupChannel object.
#' @param contFrac The contamination fraction.  Either a constant, in which case the same value is used for all cells, or a named vector, in which case the value is set for each cell.
#' @param forceAccept A warning or error is usually returned for extremely high contamination fractions.  Setting this to TRUE will turn these into messages and proceed.
#' @return A modified SoupChannel object for which the contamination (rho) has been set.
#' @examples
#' sc = load10X(system.file('extdata','toyData',package='SoupX'))
#' sc = setContaminationFraction(sc,0.1)
setContaminationFraction = function(sc,contFrac,forceAccept=FALSE){
  # Validate inputs
  validate_soup_channel(sc)
  validate_contamination_fraction(contFrac, rownames(sc$metaData))
  #Let everything else through with a diagnostic message
  if(forceAccept){
    warning = message
    stop = message
  }
  if(any(contFrac>0.5)){
    stop(sprintf("Extremely high contamination estimated (%.2g%%). This likely indicates estimation failure. ", 
                 max(contFrac)*100),
         "Check your marker genes and soup profile. Set forceAccept=TRUE to override this safety check.")
  }else if(any(contFrac>0.3)){
    warning(sprintf("Estimated contamination is very high (%.2g).",max(contFrac)))
  }
  #Now do the setting
  if(length(contFrac)==1){
    sc$metaData$rho=contFrac
  }else{
    if(!all(names(contFrac) %in% rownames(sc$metaData)))
      stop("When providing per-cell contamination fractions, 'contFrac' must be a named vector. ",
           "Names must match cell IDs (rownames of sc$metaData). ",
           "For global contamination, provide a single value. ",
           "Missing cells: ", paste(head(setdiff(names(contFrac), rownames(sc$metaData)), 10), collapse=", "))
    sc$metaData$rho[match(names(contFrac),rownames(sc$metaData))] = contFrac
  }
  return(sc)
}

#' Manually set dimension reduction for a channel
#'
#' Manually specify the dimension reduction 
#'
#' @export
#' @param sc A SoupChannel object.
#' @param DR The dimension reduction coordinates (e.g., tSNE).  This must be a data.frame, with two columns giving the two dimension reduction coordinates.  The data.frame must either have row names matching the row names of sc$metaData, or be ordered in the same order as sc$metaData.
#' @param reductName What to name the reduction (defaults to column names provided).
#' @return A modified SoupChannel object for which the dimension reduction has been set.
#' @examples
#' sc = load10X(system.file('extdata','toyData',package='SoupX'))
#' mDat = read.table(system.file('extdata','toyData','metaData.tsv',package='SoupX'),sep='\t')
#' sc = setDR(sc,mDat[,c('tSNE_1','tSNE_2')])
setDR = function(sc,DR,reductName=NULL){
  #If more than two columns, keep the first two
  if(ncol(DR)>2){
    warning(sprintf("DR has %d columns where 2 were expected.  Using first two.",ncol(DR)))
    DR = DR[,1:2]
  }
  #Check if the rownames match the metadata
  m = match(rownames(sc$metaData),rownames(DR))
  if(any(is.na(m))){
    #Can't use row-names, so the number better match
    if(nrow(DR)!=nrow(sc$metaData)){
      stop(sprintf("Rownames present in metaData not found in DR and row numbers differ (%d in metaData, %d in DR).  Each cell must have a corresponding entry in DR.",nrow(sc$metaData),nrow(DR)))
    }
    m = seq(nrow(DR))
  }
  #Should we change the names?
  if(!is.null(reductName))
    colnames(DR) = paste0(reductName,'_',1:2)
  #Add the entries in
  sc$metaData = cbind(sc$metaData,DR[m,])
  #And point them in the right place
  sc$DR = colnames(DR)
  return(sc)
}

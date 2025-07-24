#' Get expression profile of soup
#'
#' This is usually called by \code{\link{SoupChannel}},rather than directly by the user.  Uses the empty droplets in the range provided to calculate the expression profile of the soup under the assumption that these droplets only contain background.
#'
#' @export
#' @param sc A \code{SoupChannel} object.
#' @param soupRange Droplets with total UMI count in this range (excluding endpoints) are used to estimate soup.
#' @param keepDroplets Storing the full table of counts for all droplets uses a lot of space and is really only used to estimate the soup profile.  Therefore,it is dropped after the soup profile has been estimated unless this is set to \code{TRUE}.
#' @return A modified version of \code{sc} with an extra \code{soupProfile} entry containing a data.frame with the soup profile and confidence limits for all genes.
#' @examples
#' #Load droplet and count tables
#' tod <- Seurat::Read10X(system.file('extdata','toyData','raw_gene_bc_matrices','GRCh38',
#'                                   package='SoupX'))
#' toc <- Seurat::Read10X(system.file('extdata','toyData','filtered_gene_bc_matrices','GRCh38',
#'                                   package='SoupX'))
#' #Suppress calculation of soup profile automatically on load
#' sc <- SoupChannel(tod,toc,calcSoupProfile=FALSE)
#' #Retain table of droplets
#' sc <- estimateSoup(sc,keepDroplets=TRUE)
#' #Or use non-default values
#' sc <- estimateSoup(sc,soupRange=c(60,100))
estimateSoup <- function(sc,soupRange=c(0,100),keepDroplets=FALSE){
  # Validate inputs
  validate_soup_channel(sc)
  
  # Validate parameters
  if(!is.numeric(soupRange) || length(soupRange) != 2 || soupRange[1] >= soupRange[2]) {
    stop("soupRange must be a numeric vector of length 2 with increasing values. Got: [",
         paste(soupRange,collapse=","),"]")
  }
  
  # Check if soup profile already exists and tod is NULL
  if(!is.null(sc$soupProfile) && is.null(sc$tod)) {
    message("Soup profile already exists and table of droplets is NULL. Returning existing soup profile.")
    return(sc)
  }
  
  # Check if tod exists
  if(is.null(sc$tod)) {
    stop("Table of droplets (tod,call. = FALSE) is NULL. Cannot estimate soup profile without droplet data.")
  }
  
  #Estimate the soup 
  w <- which(sc$nDropUMIs > soupRange[1] & sc$nDropUMIs < soupRange[2])
  
  if(length(w) == 0) {
    stop("No droplets found in the specified soupRange [",soupRange[1],",",soupRange[2],"]. ",
         "Available UMI range: [",min(sc$nDropUMIs),",",max(sc$nDropUMIs),"]")
  }
  
  sc$soupProfile <- data.frame(row.names=rownames(sc$tod),
                              est <- rowSums(sc$tod[,w,drop=FALSE])/sum(sc$tod[,w]),
                              counts <- rowSums(sc$tod[,w,drop=FALSE]))
  
  #Saves a lot of space if we can drop the droplets now we're done with them
  if(!keepDroplets)
    sc$tod=NULL
  return(sc)
}

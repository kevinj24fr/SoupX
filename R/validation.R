#' Input validation utilities for SoupX
#'
#' Collection of validation functions to ensure input data is properly formatted
#' and contains expected values before processing.
#'
#' @name validation-utils
NULL

#' Validate SoupChannel object
#' 
#' @param sc Object to validate
#' @param require_soup_profile Logical,whether soup profile is required
#' @param require_rho Logical,whether contamination fractions are required
#' @param require_clusters Logical,whether clustering information is required
#' @return Invisibly returns TRUE if validation passes,otherwise stops with error
#' @keywords internal
validate_soup_channel <- function(sc,require_soup_profile = FALSE,require_rho = FALSE,require_clusters = FALSE) {
  if (!is(sc,'SoupChannel')) {
    stop("'sc' must be a SoupChannel object. Got object of class: ",class(sc)[1],
         ". Create SoupChannel with SoupChannel(tod,toc).",call. = FALSE)
  }
  
  if (require_soup_profile && is.null(sc$soupProfile)) {
    stop("Soup profile required but not found. Run estimateSoup(sc,call. = FALSE) first.",call. = FALSE)
  }
  
  if (require_rho && !'rho' %in% colnames(sc$metaData)) {
    stop("Contamination fractions must be calculated before this operation. ",
         "Run autoEstCont(sc) or setContaminationFraction(sc,rho) first.",call. = FALSE)
  }
  
  if (require_clusters && !'clusters' %in% colnames(sc$metaData)) {
    stop("Clustering information required for this operation. ",
         "Run setClusters(sc,cluster_vector) first,where cluster_vector assigns each cell to a cluster.",call. = FALSE)
  }
  
  invisible(TRUE)
}

#' Validate gene list for contamination estimation
#'
#' @param gene_list List or vector of genes
#' @param available_genes Character vector of available gene names  
#' @return Invisibly returns TRUE if validation passes,otherwise stops with error
#' @keywords internal
validate_gene_list <- function(gene_list,available_genes = NULL) {
  if (is.null(gene_list)) {
    stop("Gene list cannot be NULL. Provide a named list of gene sets. ",
         "Example: list(HB <- c('HBB','HBA2'),IG <- c('IGHA1','IGHG1'))",call. = FALSE)
  }
  
  if (!is.list(gene_list)) {
    stop("'gene_list' must be a named list of gene sets for contamination estimation. ",
         "Example: list(HB <- c('HBB','HBA2'),IG <- c('IGHA1','IGHG1')). ",
         "Got object of class: ",class(gene_list)[1],call. = FALSE)
  }
  
  if (length(gene_list) == 0) {
    stop("Gene list is empty. Provide at least one gene set.",call. = FALSE)
  }
  
  if (is.null(names(gene_list)) || any(names(gene_list) == "")) {
    stop("All gene sets must be named. ",
         "Example: list(HB <- c('HBB','HBA2'),IG <- c('IGHA1','IGHG1'))",call. = FALSE)
  }
  
  # Check if genes exist in available genes
  if (!is.null(available_genes)) {
    all_genes <- unique(unlist(gene_list))
    missing_genes <- setdiff(all_genes,available_genes)
    if (length(missing_genes) > 0) {
      stop("Some genes not found in data: ",
           paste(head(missing_genes,10),collapse = ","),
           if (length(missing_genes) > 10) paste(" and",length(missing_genes) - 10,"more") else "",
           call. = FALSE)
    }
  }
  
  invisible(TRUE)
}

#' Validate cluster assignments
#'
#' @param clusters Vector of cluster assignments
#' @param cell_names Character vector of expected cell names
#' @return Invisibly returns TRUE if validation passes,otherwise stops with error
#' @keywords internal
validate_clusters <- function(clusters,cell_names) {
  if (is.null(clusters)) {
    stop("Cluster assignments cannot be NULL.",call. = FALSE)
  }
  
  if (!is.character(clusters) && !is.numeric(clusters) && !is.logical(clusters)) {
    stop("Invalid cluster specification: 'clusters' must be character,numeric,or logical vector. ",
         "Got object of class: ",class(clusters)[1],call. = FALSE)
  }
  
  if (any(is.na(clusters))) {
    na_cells <- names(clusters)[is.na(clusters)]
    if (is.null(na_cells)) na_cells <- which(is.na(clusters))
    stop("NA values found in cluster assignments. All cells must be assigned to a cluster. ",
         "Cells with NA clusters: ",paste(head(na_cells,10),collapse = ","),call. = FALSE)
  }
  
  # Check that all cells are covered if clusters is named
  if (!is.null(names(clusters))) {
    missing_cells <- setdiff(cell_names,names(clusters))
    if (length(missing_cells) > 0) {
      stop("Invalid cluster specification. 'clusters' must be a named vector containing all cell names. ",
           "Missing cells: ",paste(head(missing_cells,10),collapse = ","),call. = FALSE)
    }
  } else {
    # If not named,must have same length as cell names
    if (length(clusters) != length(cell_names)) {
      stop("Invalid cluster specification: 'clusters' must be either a named vector with all cell names,",
           "or an unnamed vector with length equal to number of cells (",length(cell_names),"). ",
           "Got length: ",length(clusters),call. = FALSE)
    }
  }
  
  invisible(TRUE)
}

#' Validate contamination fractions
#'
#' @param cont_frac Numeric vector of contamination fractions
#' @param cell_names Character vector of cell names (for named vectors)
#' @return Invisibly returns TRUE if validation passes,otherwise stops with error
#' @keywords internal
validate_contamination_fraction <- function(cont_frac,cell_names = NULL) {
  if (!is.numeric(cont_frac)) {
    stop("Contamination fractions must be numeric. Got class: ",class(cont_frac)[1],call. = FALSE)
  }
  
  if (any(is.na(cont_frac))) {
    stop("Contamination fractions cannot contain NA values.",call. = FALSE)
  }
  
  if (any(cont_frac < 0)) {
    stop("Contamination fractions must be non-negative. ",
         "Minimum value: ",min(cont_frac),call. = FALSE)
  }
  
  if (any(cont_frac > 1)) {
    stop("Contamination fraction greater than 1.0 detected (impossible,call. = FALSE). ",
         "Maximum value: ",max(cont_frac),". This indicates an error in estimation. ",
         "Check your soup profile and marker gene selection.",call. = FALSE)
  }
  
  # If providing per-cell fractions,must be named
  if (length(cont_frac) > 1 && !is.null(cell_names)) {
    if (is.null(names(cont_frac))) {
      stop("When providing per-cell contamination fractions,'cont_frac' must be a named vector. ",
           "Names must match cell IDs. For global contamination,provide a single value.",call. = FALSE)
    }
    
    missing_cells <- setdiff(names(cont_frac),cell_names)
    if (length(missing_cells) > 0) {
      stop("Cell names in contamination fractions don't match expected cells. ",
           "Missing cells: ",paste(head(missing_cells,10),collapse = ","),call. = FALSE)
    }
  }
  
  invisible(TRUE)
}

#' Validate dimension reduction matrix
#'
#' @param DR Matrix or data.frame with dimension reduction coordinates
#' @param cell_names Character vector of expected cell names
#' @param min_dims Minimum number of dimensions required
#' @return Invisibly returns TRUE if validation passes,otherwise stops with error
#' @keywords internal
validate_dimension_reduction <- function(DR,cell_names,min_dims = 2) {
  if (is.null(DR)) {
    stop("Dimension reduction matrix cannot be NULL.",call. = FALSE)
  }
  
  DR <- as.data.frame(DR)
  
  if (ncol(DR) < min_dims) {
    stop("Dimension reduction matrix 'DR' must have at least ",min_dims," columns for plotting. ",
         "Got ",ncol(DR)," columns.",call. = FALSE)
  }
  
  if (nrow(DR) == 0) {
    stop("Dimension reduction matrix is empty.",call. = FALSE)
  }
  
  if (!is.null(cell_names)) {
    missing_cells <- setdiff(rownames(DR),cell_names)
    if (length(missing_cells) > 0) {
      stop("Rownames of 'DR' must match cell names in count matrix. ",
           "Missing cells: ",paste(head(missing_cells,10),collapse = ","),call. = FALSE)
    }
  }
  
  invisible(TRUE)
}

#' Validate count matrices for compatibility
#'
#' @param tod Table of droplets (full count matrix)
#' @param toc Table of counts (cell-only matrix)
#' @return Invisibly returns TRUE if validation passes,otherwise stops with error
#' @keywords internal
validate_count_matrices <- function(tod,toc) {
  if (is.null(tod) || is.null(toc)) {
    stop("Count matrices cannot be NULL.",call. = FALSE)
  }
  
  if (nrow(tod) != nrow(toc)) {
    stop("Gene mismatch: table of droplets (tod,call. = FALSE) has ",nrow(tod)," genes but table of counts (toc) has ",nrow(toc)," genes. ",
         "Both matrices must have the same genes in the same order.",call. = FALSE)
  }
  
  if (any(rownames(tod) != rownames(toc))) {
    first_diff <- which(rownames(tod) != rownames(toc))[1]
    stop("Gene order mismatch: rownames of tod and toc differ. ",
         "Both matrices must have identical gene names in the same order. ",
         "First difference at position: ",first_diff,call. = FALSE)
  }
  
  if (ncol(toc) > ncol(tod)) {
    stop("Table of counts (toc,call. = FALSE) cannot have more cells than table of droplets (tod). ",
         "toc: ",ncol(toc)," cells,tod: ",ncol(tod)," droplets.",call. = FALSE)
  }
  
  invisible(TRUE)
}

#' Validate soup profile
#'
#' @param soup_profile Data.frame with soup profile information
#' @param gene_names Character vector of expected gene names
#' @return Invisibly returns TRUE if validation passes,otherwise stops with error
#' @keywords internal
validate_soup_profile <- function(soup_profile,gene_names = NULL) {
  if (is.null(soup_profile)) {
    stop("Soup profile cannot be NULL.",call. = FALSE)
  }
  
  if (!is.data.frame(soup_profile)) {
    stop("Soup profile must be a data.frame. Got class: ",class(soup_profile)[1],call. = FALSE)
  }
  
  required_cols <- c("est","counts")
  missing_cols <- setdiff(required_cols,colnames(soup_profile))
  if (length(missing_cols) > 0) {
    stop("Soup profile missing required columns: ",paste(missing_cols,collapse = ",",call. = FALSE),". ",
         "Required columns: ",paste(required_cols,collapse = ","),". ",
         "Current columns: ",paste(colnames(soup_profile),collapse = ","),call. = FALSE)
  }
  
  if (!is.null(gene_names)) {
    missing_genes <- setdiff(gene_names,rownames(soup_profile))
    if (length(missing_genes) > 0) {
      stop("Soup profile missing genes found in count matrix. ",
           "Missing genes: ",paste(head(missing_genes,10),collapse = ","),call. = FALSE)
    }
  }
  
  # Validate data ranges
  if (any(soup_profile$est < 0,na.rm = TRUE)) {
    stop("Soup profile estimates cannot be negative.",call. = FALSE)
  }
  
  if (any(soup_profile$counts < 0,na.rm = TRUE)) {
    stop("Soup profile counts cannot be negative.",call. = FALSE)
  }
  
  invisible(TRUE)
}

#' Validate numeric parameters for SoupX functions
#'
#' @param value Numeric value to validate
#' @param param_name Name of parameter for error messages
#' @param min_value Minimum allowed value (inclusive)
#' @param max_value Maximum allowed value (inclusive)
#' @param allow_na Whether NA values are allowed
#' @return Invisibly returns TRUE if validation passes,otherwise stops with error
#' @keywords internal
validate_numeric_parameter <- function(value,param_name,min_value = -Inf,max_value = Inf,allow_na = FALSE) {
  if (!allow_na && any(is.na(value))) {
    stop("Parameter '",param_name,"' cannot contain NA values.",call. = FALSE)
  }
  
  if (!is.numeric(value)) {
    stop("Parameter '",param_name,"' must be numeric. Got class: ",class(value)[1],call. = FALSE)
  }
  
  if (any(value[!is.na(value)] < min_value,na.rm = TRUE)) {
    stop("Parameter '",param_name,"' must be >= ",min_value,". ",
         "Minimum value found: ",min(value,na.rm = TRUE),call. = FALSE)
  }
  
  if (any(value[!is.na(value)] > max_value,na.rm = TRUE)) {
    stop("Parameter '",param_name,"' must be <= ",max_value,". ",
         "Maximum value found: ",max(value,na.rm = TRUE),call. = FALSE)
  }
  
  invisible(TRUE)
}
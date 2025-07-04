\name{validation-utils}
\alias{validation-utils}
\alias{validate_soup_channel}
\alias{validate_gene_list}
\alias{validate_clusters}
\alias{validate_contamination_fraction}
\alias{validate_dimension_reduction}
\alias{validate_count_matrices}
\alias{validate_soup_profile}
\title{Input Validation Utilities for SoupX}
\description{
Collection of validation functions to ensure input data is properly formatted 
and contains expected values before processing. These are internal functions 
used throughout the package to provide consistent error checking and helpful 
error messages.
}
\details{
These validation functions are used internally to:
\itemize{
  \item Check object types and structure
  \item Validate parameter ranges and formats
  \item Ensure data consistency between related objects
  \item Provide informative error messages for common user mistakes
}

The main validation functions include:
\itemize{
  \item \code{validate_soup_channel}: Validates SoupChannel objects
  \item \code{validate_gene_list}: Validates gene lists for contamination estimation
  \item \code{validate_clusters}: Validates cluster assignments
  \item \code{validate_contamination_fraction}: Validates contamination fractions
  \item \code{validate_dimension_reduction}: Validates dimension reduction matrices
  \item \code{validate_count_matrices}: Validates count matrix compatibility
  \item \code{validate_soup_profile}: Validates soup profile data frames
}
}
\keyword{internal} 
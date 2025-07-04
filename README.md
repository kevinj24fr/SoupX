# SoupX - Enhanced Version

An enhanced R package for the estimation and removal of cell free mRNA contamination in droplet based single cell RNA-seq data.

This fork includes significant improvements over the original repository, focusing on stability, reliability, and developer experience.

## 🚀 Major Improvements

### ✅ **Enhanced Stability & Reliability**
- **Comprehensive validation framework** - Input validation for all functions
- **Improved error handling** - Clear, informative error messages
- **Better input validation** - Prevents crashes from invalid inputs
- **Robust matrix handling** - Fixed issues with unnamed matrices

### ✅ **Testing & Quality Assurance**
- **Complete test suite** - Comprehensive testing with testthat
- **Real data testing** - Includes PBMC3k dataset for validation
- **Test scripts** - Ready-to-use testing scripts for validation
- **Quality metrics** - Performance and reliability tracking

### ✅ **Developer Experience**
- **Modern dependencies** - Updated package dependencies
- **Better documentation** - Enhanced function documentation
- **Debugging tools** - Diagnostic scripts for troubleshooting
- **Installation improvements** - Streamlined setup process

### ✅ **Key Bug Fixes**
- **Matrix naming fix** - Resolved issues with 10x data loading when matrices lack proper row/column names
- **Ensembl ID support** - Proper handling of gene identifiers
- **Memory optimization** - Improved memory usage for large datasets

## Installation

Install the enhanced version with all improvements:

```R
devtools::install_github("kevinj24fr/SoupX")
```

If you encounter errors saying `multtest` is unavailable, please install this manually from bioconductor with:

```R
BiocManager::install('multtest')
```

## Quickstart

Decontaminate one channel of 10X data mapped with cellranger by running:

```R
sc = load10X('path/to/your/cellranger/outs/folder')
sc = autoEstCont(sc)
out = adjustCounts(sc)
```

or to manually load decontaminate any other data

```R
sc = SoupChannel(table_of_droplets,table_of_counts)
sc = setClusters(sc,cluster_labels)
sc = autoEstCont(sc)
out = adjustCounts(sc)
```

`out` will then contain a corrected matrix to be used in place of the original table of counts in downstream analyses.

## Testing the Package

This fork includes comprehensive testing capabilities:

```R
# Run the test suite
library(testthat)
test_package("SoupX")

# Test with real 10x data (included in package)
source(system.file("test_data", "debug_soupx_fixed.R", package = "SoupX"))
```

## Important Notes for 10x Data

**⚠️ Critical Fix:** When loading 10x data manually, ensure your matrices have proper row and column names:

```R
# Load data
raw_data <- readMM('raw_gene_bc_matrices/hg19/matrix.mtx')
filtered_data <- readMM('filtered_gene_bc_matrices/hg19/matrix.mtx')
genes <- read.table('raw_gene_bc_matrices/hg19/genes.tsv', header=FALSE, stringsAsFactors=FALSE)
barcodes_raw <- read.table('raw_gene_bc_matrices/hg19/barcodes.tsv', header=FALSE, stringsAsFactors=FALSE)
barcodes_filtered <- read.table('filtered_gene_bc_matrices/hg19/barcodes.tsv', header=FALSE, stringsAsFactors=FALSE)

# Set proper names (use Ensembl IDs for row names to avoid duplicates)
rownames(raw_data) <- genes$V1  # Ensembl IDs
rownames(filtered_data) <- genes$V1
colnames(raw_data) <- barcodes_raw$V1
colnames(filtered_data) <- barcodes_filtered$V1

# Create SoupChannel
sc <- SoupChannel(raw_data, filtered_data, calcSoupProfile=TRUE)
```

## Documentation

The methodology implemented in this package is explained in detail in [this paper](https://doi.org/10.1093/gigascience/giaa151).  

A detailed vignette is provided with the package and can be viewed [here](https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html).  

## Citing SoupX

If you use SoupX in your work, please cite: "Young, M.D., Behjati, S. (2020). SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data, GigaScience, Volume 9, Issue 12, December 2020, giaa151bioRxiv, 303727, https://doi.org/10.1093/gigascience/giaa151"

## Changelog

### v1.6.1 (This Enhanced Fork)

**Major Improvements:**
- ✅ Added comprehensive validation framework (`R/validation.R`)
- ✅ Implemented robust error handling throughout the package
- ✅ Created complete test suite with testthat
- ✅ Fixed matrix naming issues for 10x data loading
- ✅ Added real data testing with PBMC3k dataset
- ✅ Updated dependencies and package metadata
- ✅ Enhanced function documentation
- ✅ Added debugging and diagnostic tools
- ✅ Improved memory efficiency for large datasets

**Bug Fixes:**
- Fixed issues with unnamed matrices causing soup profile calculation failures
- Resolved duplicate row name errors when using gene symbols
- Improved handling of Ensembl IDs and gene identifiers
- Enhanced input validation to prevent crashes

**Code Quality:**
- Removed commented out and unnecessary code
- Cleaned up legacy documentation
- Streamlined package structure
- Improved code organization and readability

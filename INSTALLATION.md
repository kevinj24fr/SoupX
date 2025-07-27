# SoupX Installation Guide

## Quick Installation

### For Academic Use (Free)

```r
# Install from GitHub
devtools::install_github('kevinj24fr/SoupX')

# Load the package
library(SoupX)
```

### For Commercial Use (Requires License)

Contact: kevinj24fr@gmail.com

## System Requirements

- R version 4.0 or higher
- Matrix package (>= 1.3.0)
- Optional: Seurat (>= 4.0.0) for 10X data loading
- Optional: ggplot2 (>= 3.0.0) for visualization

## Dependencies

### Required Dependencies
- **Matrix**: For sparse matrix operations
- **methods**: For S3 class methods

### Suggested Dependencies
- **Seurat**: For loading 10X Genomics data
- **ggplot2**: For plotting functions
- **knitr**: For vignette building
- **testthat**: For running tests

## Installation Troubleshooting

### If you encounter compilation errors:

1. **Update R and packages**:
```r
update.packages()
```

2. **Install from source if needed**:
```r
# Install devtools first
install.packages("devtools")

# Then install SoupX
devtools::install_github('kevinj24fr/SoupX', dependencies = TRUE)
```

3. **For macOS users with compilation issues**:
```r
# Install Xcode command line tools
# xcode-select --install

# Then install with specific repository
devtools::install_github('kevinj24fr/SoupX', repos = "https://cran.r-project.org")
```

## Verification

After installation, verify the package works:

```r
library(SoupX)

# Load test data
data(scToy)

# Test basic functionality
sc <- scToy
sc <- estimateSoup(sc)
sc <- setContaminationFraction(sc, 0.1)
adjusted <- adjustCounts(sc)

cat("Installation successful! Package is ready to use.\n")
```

## Usage Examples

### Basic Usage
```r
# Load your data
sc <- load10X("path/to/your/10x/data")

# Estimate contamination
sc <- autoEstCont(sc)

# Adjust counts
adjusted <- adjustCounts(sc)
```

### For PBMC Data
See `pbmc_soupx_analysis.R` for comprehensive PBMC analysis examples.

## Support

- **GitHub Issues**: https://github.com/kevinj24fr/SoupX/issues
- **Documentation**: See `FINAL_TEST_SUMMARY.md` for detailed test results
- **Commercial Licensing**: Contact kevinj24fr@gmail.com

## Version History

- **v1.6.6**: Production ready with comprehensive testing
- **v1.6.5**: Critical bug fixes and performance optimizations
- **v1.6.4**: Enhanced error handling and validation

---

**Ready for Production Use!** 🚀 
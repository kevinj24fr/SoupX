# SoupX - Robust Single Cell RNA-seq Contamination Correction

> **DUAL LICENSE:** Free for academic use • Commercial license required for business use  
> **Commercial licensing:** kevin.joseph@uniklinik-freiburg.de

---

## Release Notes: Version 1.6.4

- **FIXED:** All assignment and argument errors (`<-` vs `=`) in function calls, data.frame, and list construction.
- **FIXED:** All subsetting errors (`drop <- FALSE` vs `drop = FALSE`).
- **FIXED:** Validation logic for all numeric parameters and input checks.
- **FIXED:** Documentation and examples updated for clarity and accuracy.
- **IMPROVED:** Optional dependencies (Seurat, ggplot2) are handled gracefully.
- **ENHANCED:** Comprehensive error handling with actionable messages.
- **ADDED:** Edge case and performance regression tests.
- **STANDARDIZED:** Code style and documentation across all functions.
- **OPTIMIZED:** Memory usage and sparse matrix operations.
- **EXPANDED:** Test coverage for all core and edge-case scenarios.

---

## Licensing

**SoupX uses a dual licensing model:**

### **Academic Use (FREE)**
- Universities and research institutions
- Educational use in academic courses  
- Non-profit research organizations
- Personal academic research
- Open source academic projects

### **Commercial Use (LICENSE REQUIRED)**
- Biotechnology companies
- Pharmaceutical companies  
- Commercial research organizations
- Consulting services
- SaaS platforms
- Any for-profit use

---

## Installation

### Academic Users (Free)
```r
# Install from GitHub (academic use only)
if (!require(devtools)) install.packages("devtools")
devtools::install_github("kevinj24fr/SoupX")

# Or install from local source
devtools::install_local("path/to/SoupX")

# Optional dependencies for enhanced functionality
install.packages(c("Seurat", "ggplot2"))  # For 10X loading and plotting
```

## Key Features

### Flexible Dependencies
```r
# Works without optional packages
library(SoupX)
sc <- SoupChannel(tod, toc)  # Core functionality always available

# Enhanced features with optional packages
if (requireNamespace("Seurat", quietly = TRUE)) {
  sc <- load10X("path/to/10x/data")  # 10X data loading
}

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plotSoupCorrelation(sc)  # Visualization functions
}
```

### Robust Error Handling
```r
# Clear, actionable error messages
sc <- setContaminationFraction(sc, 1.5)
# Error: Contamination fraction greater than 1.0 detected (impossible).
# Maximum value: 1.5. This indicates an error in estimation.
# Check your soup profile and marker gene selection.

# Graceful handling of edge cases
sc <- autoEstCont(homogeneous_data)
# Error: No suitable marker genes found for contamination estimation.
# Try: (1) reducing tfidfMin to accept less specific markers,
# or (2) reducing soupQuantile to include lower-expressed genes,
# or (3) manually set contamination with setContaminationFraction().
```

## 📖 Quick Start

```r
library(SoupX)

# Load 10x data
sc <- load10X("path/to/10x/data")

# Estimate contamination
sc <- autoEstCont(sc)

# Correct counts
adjusted <- adjustCounts(sc)

# Generate comprehensive report
report <- generateSoupXReport(sc, adjusted)
```

## Performance Benchmarking

```r
# Benchmark performance
benchmark_results <- benchmark_soupx(sc, iterations = 5)
print(benchmark_results)
```

## Testing

```r
# Run test suite
devtools::test()

# Run comprehensive test with real data
source("test_data/visualization_examples.R")
```

## Visualization Examples

See `test_data/visualization_examples.R` for comprehensive examples of all new visualization capabilities.

## Changelog

### Version 1.6.4 (Current) - Robustness & Quality
- **FIXED**: Critical syntax error in core functions  
- **IMPROVED**: Optional dependencies (Seurat, ggplot2) with graceful fallbacks
- **ENHANCED**: Comprehensive error handling with actionable messages
- **ADDED**: Edge case testing and validation framework
- **STANDARDIZED**: Code style and documentation across all functions
- **OPTIMIZED**: Memory usage and sparse matrix operations
- **EXPANDED**: Test coverage including performance regression tests

### Version 1.6.3 
- **NEW**: Comprehensive visualization suite
- **NEW**: Quality control dashboard and reporting
- **IMPROVED**: Performance optimizations and benchmarking tools

### Version 1.6.2
- **ADDED**: Input validation framework and enhanced error messages
- **IMPROVED**: Memory-efficient sparse matrix operations

### Version 1.6.1 
- **FIXED**: Backward compatibility and 10x data loading issues
- **ADDED**: Comprehensive test suite

## Acknowledgments

Original SoupX by Matthew Young. Enhanced with robustness improvements, improved error handling, comprehensive validation, and production-ready architecture for reliable single-cell RNA-seq contamination correction.

---

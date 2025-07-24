# SoupX v1.6.4 - Robust, Bug-Free Release

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

## Key Improvements

This release represents a comprehensive overhaul of the SoupX package:

- **Robust Error Handling:** All functions now provide clear, actionable error messages
- **Dual Licensing:** Free for academic use, commercial license required for business use
- **Comprehensive Testing:** Edge cases and performance regression tests included
- **Production Ready:** All known bugs fixed, syntax errors corrected
- **Enhanced Documentation:** Updated examples and clear usage instructions

## Installation

Academic users can install directly from GitHub:
```r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("kevinj24fr/SoupX")
```

Commercial users require a license. Contact: kevin.joseph@uniklinik-freiburg.de 
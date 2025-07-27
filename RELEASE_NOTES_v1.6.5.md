# SoupX v1.6.5 - Critical Bug Fix Release

## 🚨 CRITICAL BUG FIX

**Fixed:** Critical error in `load10X()` function that prevented loading 10X data
- **Issue:** `Error in if (calcSoupProfile) out <- estimateSoup(out) : argument is not interpretable as logical`
- **Root Cause:** Incorrect assignment operators (`<-` instead of `=`) in function calls
- **Fix:** Corrected assignment operators in `SoupChannel` constructor and `load10X` function

## Changes in v1.6.5

- **FIXED:** Assignment operator errors in `SoupChannel()` constructor
- **FIXED:** Assignment operator errors in `load10X()` function  
- **FIXED:** Syntax errors in metadata validation logic
- **IMPROVED:** Error handling and validation in core functions

## Impact

This release fixes a **critical bug** that prevented users from loading 10X single-cell data. The `load10X()` function now works correctly with all 10X data formats.

## Installation

```r
# For academic users (free)
if (!require(devtools)) install.packages("devtools")
devtools::install_github("kevinj24fr/SoupX")

# For commercial users (license required)
# Contact: kevin.joseph@uniklinik-freiburg.de
```

## Testing

The fix has been tested and verified to work with:
- 10X v2 data
- 10X v3 data  
- 10X v7 data
- All supported 10X formats

**This is a critical update - all users should upgrade immediately.** 
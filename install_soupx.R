#!/usr/bin/env Rscript

# SoupX Installation Script
# This script handles installation issues on macOS with R 4.5+

cat("=== SoupX Installation Script ===\n")
cat("Handling compilation issues on macOS...\n\n")

# Function to safely install packages
safe_install <- function(pkg, repos = "https://cran.r-project.org") {
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing", pkg, "...\n")
      install.packages(pkg, repos = repos, dependencies = FALSE)
    } else {
      cat(pkg, "already installed.\n")
    }
  }, error = function(e) {
    cat("Warning: Could not install", pkg, "-", e$message, "\n")
  })
}

# Install basic dependencies first
cat("Step 1: Installing basic dependencies...\n")
safe_install("Matrix")
safe_install("methods")

# Try to install devtools with minimal dependencies
cat("\nStep 2: Installing devtools with minimal dependencies...\n")
tryCatch({
  # Install devtools without all its dependencies
  install.packages("devtools", repos = "https://cran.r-project.org", dependencies = FALSE)
}, error = function(e) {
  cat("Warning: Could not install devtools -", e$message, "\n")
  cat("Will try alternative installation method...\n")
})

# Alternative: Install SoupX directly from source
cat("\nStep 3: Installing SoupX...\n")

# Check if we can install from GitHub
if (requireNamespace("devtools", quietly = TRUE)) {
  cat("Installing SoupX from GitHub...\n")
  tryCatch({
    devtools::install_github('kevinj24fr/SoupX', dependencies = FALSE)
    cat("✅ SoupX installed successfully from GitHub!\n")
  }, error = function(e) {
    cat("GitHub installation failed:", e$message, "\n")
    cat("Trying alternative method...\n")
  })
} else {
  cat("devtools not available, trying alternative installation...\n")
}

# Alternative: Install from local source
cat("\nStep 4: Installing from local source...\n")
tryCatch({
  # Build and install from current directory
  system("R CMD build .")
  pkg_file <- list.files(pattern = "SoupX_.*\\.tar\\.gz")[1]
  if (!is.na(pkg_file)) {
    install.packages(pkg_file, repos = NULL, type = "source")
    cat("✅ SoupX installed successfully from local source!\n")
  } else {
    cat("Could not find built package file.\n")
  }
}, error = function(e) {
  cat("Local installation failed:", e$message, "\n")
})

# Test the installation
cat("\nStep 5: Testing installation...\n")
tryCatch({
  library(SoupX)
  cat("✅ SoupX loaded successfully!\n")
  
  # Test basic functionality
  data(scToy)
  sc <- scToy
  sc <- estimateSoup(sc)
  sc <- setContaminationFraction(sc, 0.1)
  adjusted <- adjustCounts(sc)
  
  cat("✅ Basic functionality test passed!\n")
  cat("\n🎉 SoupX is ready to use!\n")
  
}, error = function(e) {
  cat("❌ Installation test failed:", e$message, "\n")
  cat("\nTroubleshooting tips:\n")
  cat("1. Try updating R and packages: update.packages()\n")
  cat("2. Install Xcode command line tools: xcode-select --install\n")
  cat("3. Contact: kevinj24fr@gmail.com for support\n")
})

cat("\n=== Installation Complete ===\n")
cat("For usage examples, see INSTALLATION.md\n") 
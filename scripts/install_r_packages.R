#!/usr/bin/env Rscript

####################################################################################
# Script: install_r_packages.R
# Purpose: Install all required R packages for IISAGE ATAC-seq analysis
# Author: Louis Paul Decena-Segarra, PhD
# Usage: Rscript install_r_packages.R
#        Or from R: source("install_r_packages.R")
####################################################################################

cat("\n")
cat("========================================================================\n")
cat("  IISAGE ATAC-seq R Package Installation\n")
cat("========================================================================\n\n")

# Clean up any stale lock files
cat("Cleaning up any stale lock files...\n")
lib_path <- .libPaths()[1]
lock_files <- list.files(lib_path, pattern = "^00LOCK-", full.names = TRUE)
if (length(lock_files) > 0) {
  cat("Removing", length(lock_files), "lock files:\n")
  for (lock in lock_files) {
    cat("  -", basename(lock), "\n")
    unlink(lock, recursive = TRUE)
  }
} else {
  cat("No lock files found.\n")
}
cat("\n")

# Track installation status
install_log <- list(
  success = character(),
  failed = character(),
  skipped = character()
)

# ——————————————————————————————————————————————————————————————————————————————
# Function to install and test a package
# ——————————————————————————————————————————————————————————————————————————————

install_and_test <- function(pkg, source = "CRAN", required = TRUE) {
  pkg_label <- ifelse(required, "[REQUIRED]", "[OPTIONAL]")
  cat(sprintf("%-50s ... ", paste0(pkg, " ", pkg_label)))

  # Check if already installed
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✓ Already installed\n")
    install_log$skipped <<- c(install_log$skipped, pkg)
    return(TRUE)
  }

  # Try to install
  tryCatch({
    if (source == "CRAN") {
      install.packages(pkg, repos = "https://cloud.r-project.org",
                      dependencies = TRUE, quiet = TRUE)
    } else if (source == "Bioconductor") {
      BiocManager::install(pkg, update = FALSE, ask = FALSE, force = TRUE,
                          type = "source", quiet = FALSE)
    }

    # Test if installation worked
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat("✓ Installed successfully\n")
      install_log$success <<- c(install_log$success, pkg)
      return(TRUE)
    } else {
      cat("✗ Installation failed (package not loadable)\n")
      install_log$failed <<- c(install_log$failed, pkg)
      return(FALSE)
    }
  }, error = function(e) {
    cat("✗ Installation error:", conditionMessage(e), "\n")
    install_log$failed <<- c(install_log$failed, pkg)
    return(FALSE)
  })
}

# ——————————————————————————————————————————————————————————————————————————————
# Step 1: Install BiocManager (required for Bioconductor packages)
# ——————————————————————————————————————————————————————————————————————————————

cat("\n[Step 1/4] Installing BiocManager\n")
cat("========================================================================\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager ... ")
  install.packages("BiocManager", repos = "https://cloud.r-project.org",
                  quiet = TRUE)
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    cat("✓ Success\n")
  } else {
    cat("✗ Failed\n")
    stop("BiocManager installation failed. Cannot proceed.")
  }
} else {
  cat("BiocManager ... ✓ Already installed\n")
}

cat("BiocManager version:", as.character(BiocManager::version()), "\n")

# ——————————————————————————————————————————————————————————————————————————————
# Step 2: Install Core CRAN Packages (Required)
# ——————————————————————————————————————————————————————————————————————————————

cat("\n[Step 2/4] Installing Core CRAN Packages\n")
cat("========================================================================\n")

cran_packages <- c(
  "pheatmap",      # Enhanced heatmaps
  "readr",         # Fast CSV reading
  "dplyr",         # Data manipulation
  "ggplot2",       # Plotting
  "RColorBrewer",  # Color palettes
  "gridExtra"      # Multiple plots
)

for (pkg in cran_packages) {
  install_and_test(pkg, source = "CRAN", required = TRUE)
}

# ——————————————————————————————————————————————————————————————————————————————
# Step 3: Install Core Bioconductor Packages (Required)
# ——————————————————————————————————————————————————————————————————————————————

cat("\n[Step 3/4] Installing Core Bioconductor Packages\n")
cat("========================================================================\n")

bioc_core_packages <- c(
  "DiffBind",           # Differential binding analysis
  "GenomicRanges",      # Genomic interval operations
  "GenomicFeatures",    # Genomic feature handling
  "rtracklayer",        # Import/export genomic data
  "Rsamtools",          # BAM file handling
  "GenomeInfoDb"        # Genome information
)

for (pkg in bioc_core_packages) {
  install_and_test(pkg, source = "Bioconductor", required = TRUE)
}

# ——————————————————————————————————————————————————————————————————————————————
# Step 4: Install Annotation/Analysis Packages (Highly Recommended)
# ——————————————————————————————————————————————————————————————————————————————

cat("\n[Step 4/4] Installing Annotation & Analysis Packages\n")
cat("========================================================================\n")

bioc_annotation_packages <- c(
  "ChIPseeker",         # Peak annotation
  "clusterProfiler",    # GO/pathway enrichment (bonus)
  "AnnotationDbi",      # Annotation database interface
  "biomaRt"            # BioMart interface for gene info
)

for (pkg in bioc_annotation_packages) {
  install_and_test(pkg, source = "Bioconductor", required = FALSE)
}

# ——————————————————————————————————————————————————————————————————————————————
# Step 5: Install Optional Packages (For advanced features)
# ——————————————————————————————————————————————————————————————————————————————

cat("\n[Optional] Installing Advanced Analysis Packages\n")
cat("========================================================================\n")
cat("(These are optional - script will work without them)\n\n")

optional_packages <- c(
  "memes",             # Motif analysis (requires MEME suite)
  "TFBSTools",         # Transcription factor binding sites
  "motifmatchr",       # Motif matching
  "chromVAR",          # Chromatin variability
  "ATACseqQC",         # ATAC-seq specific QC
  "GenomicAlignments", # Required for ATACseqQC
  "ggseqlogo"          # Sequence logo plots
)

for (pkg in optional_packages) {
  install_and_test(pkg, source = "Bioconductor", required = FALSE)
}

# ——————————————————————————————————————————————————————————————————————————————
# Summary Report
# ——————————————————————————————————————————————————————————————————————————————

cat("\n")
cat("========================================================================\n")
cat("  Installation Summary\n")
cat("========================================================================\n\n")

cat("✓ Successfully installed:", length(install_log$success), "packages\n")
if (length(install_log$success) > 0) {
  cat("  -", paste(install_log$success, collapse = "\n  - "), "\n\n")
}

cat("⊘ Already installed (skipped):", length(install_log$skipped), "packages\n")
if (length(install_log$skipped) > 0) {
  cat("  -", paste(install_log$skipped, collapse = "\n  - "), "\n\n")
}

if (length(install_log$failed) > 0) {
  cat("✗ Failed to install:", length(install_log$failed), "packages\n")
  cat("  -", paste(install_log$failed, collapse = "\n  - "), "\n\n")
  cat("Note: Some failures are expected for optional packages.\n")
  cat("      Core functionality will work without them.\n\n")
}

# ——————————————————————————————————————————————————————————————————————————————
# Test Core Functionality
# ——————————————————————————————————————————————————————————————————————————————

cat("========================================================================\n")
cat("  Testing Core Package Loading\n")
cat("========================================================================\n\n")

test_packages <- c("DiffBind", "pheatmap", "readr", "dplyr", "ggplot2",
                   "GenomicRanges", "GenomicFeatures")

all_loaded <- TRUE
for (pkg in test_packages) {
  cat(sprintf("%-30s ... ", pkg))
  test_result <- tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    cat("✓ OK\n")
    TRUE
  }, error = function(e) {
    cat("✗ FAILED:", conditionMessage(e), "\n")
    FALSE
  })

  if (!test_result) all_loaded <- FALSE
}

cat("\n")

# ——————————————————————————————————————————————————————————————————————————————
# Final Status
# ——————————————————————————————————————————————————————————————————————————————

if (all_loaded) {
  cat("========================================================================\n")
  cat("  ✓ SUCCESS! All core packages installed and tested.\n")
  cat("========================================================================\n\n")
  cat("You can now run:\n")
  cat("  Rscript generate_atac_figs_v2.R --base_dir=/path/to/projects\n\n")

  # Check for optional enhancements
  if (requireNamespace("ChIPseeker", quietly = TRUE)) {
    cat("✓ ChIPseeker available - Peak annotation enabled\n")
  } else {
    cat("⚠ ChIPseeker not available - Peak annotation will be skipped\n")
    cat("  Install with: BiocManager::install('ChIPseeker')\n")
  }

  if (requireNamespace("memes", quietly = TRUE)) {
    cat("✓ memes available - R-based motif analysis enabled\n")
  } else {
    cat("⚠ memes not available - R-based motif analysis disabled\n")
    cat("  (HOMER-based motif analysis will still work if HOMER is installed)\n")
  }

  cat("\nFor HOMER motif analysis (optional but recommended):\n")
  cat("  Install from: http://homer.ucsd.edu/homer/\n\n")

} else {
  cat("========================================================================\n")
  cat("  ✗ WARNING: Some core packages failed to load.\n")
  cat("========================================================================\n\n")
  cat("Please check the error messages above and try:\n")
  cat("  1. Update R to the latest version\n")
  cat("  2. Update Bioconductor: BiocManager::install()\n")
  cat("  3. Reinstall failed packages manually\n\n")
}

# ——————————————————————————————————————————————————————————————————————————————
# System Information
# ——————————————————————————————————————————————————————————————————————————————

cat("========================================================================\n")
cat("  System Information\n")
cat("========================================================================\n")
cat("R version:", as.character(getRversion()), "\n")
cat("Bioconductor version:", as.character(BiocManager::version()), "\n")
cat("Platform:", R.version$platform, "\n")
cat("OS:", Sys.info()["sysname"], Sys.info()["release"], "\n")

# Save session info
session_file <- "package_installation_session_info.txt"
cat("\nSaving detailed session info to:", session_file, "\n")
sink(session_file)
cat("IISAGE ATAC-seq Package Installation\n")
cat("Date:", as.character(Sys.time()), "\n\n")
sessionInfo()
sink()

cat("\n========================================================================\n")
cat("  Installation script completed!\n")
cat("========================================================================\n\n")

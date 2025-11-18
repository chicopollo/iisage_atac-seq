#!/usr/bin/env Rscript

####################################################################################
# Script: generate_atac_figs_v2.R
# Purpose: Enhanced ATAC-seq downstream analysis with motif discovery and annotation
# Author: Louis Paul Decena-Segarra, PhD
# Version: 2.0.0
# Created: 2025-11-14
#
# New features in v2.0:
#   - Motif analysis (HOMER integration + R-based analysis)
#   - Genomic feature annotation (ChIPseeker)
#   - Comprehensive peak feature tables
#   - Distance to TSS analysis
#   - Gene annotation
#   - Summary statistics
#   - HTML reports
#
# Requirements:
#   - DiffBind, ChIPseeker, GenomicFeatures, TxDb.* (organism-specific)
#   - memes (optional, for motif analysis)
#   - HOMER (optional, system-level motif discovery)
#   - pheatmap, readr, dplyr, ggplot2, plotly
#
# Usage:
#   Rscript generate_atac_figs_v2.R --base_dir=/path/to/atac_projects
#   Or from R:
#   source("generate_atac_figs_v2.R")
#   generate_comprehensive_analysis("/path/to/atac_projects")
####################################################################################
 

  # ============================================================================
  # Minimal Auto-Install 
  # ============================================================================

  # Install if missing function
  install_if_missing <- function(pkg, source = "CRAN") {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      if (source == "CRAN") {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      } else {
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      }
    }
  }

  # Core CRAN packages (required)
  sapply(c("pheatmap", "readr", "dplyr", "ggplot2", "RColorBrewer", "gridExtra"),
         install_if_missing, source = "CRAN")

  # Core Bioconductor packages (required)
  sapply(c("DiffBind", "GenomicRanges", "GenomicFeatures", "rtracklayer",
           "Rsamtools", "GenomeInfoDb"),
         install_if_missing, source = "Bioc")

  # Annotation packages (recommended)
  sapply(c("ChIPseeker", "clusterProfiler", "AnnotationDbi", "biomaRt"),
         function(x) tryCatch(install_if_missing(x, "Bioc"), error = function(e) NULL))

  # Optional packages (advanced features)
  sapply(c("memes", "TFBSTools", "motifmatchr", "chromVAR", "ATACseqQC",
           "GenomicAlignments", "ggseqlogo"),
         function(x) tryCatch(install_if_missing(x, "Bioc"), error = function(e) NULL))

  # Clean up
  rm(install_if_missing)

  # ============================================================================


  # ——————————————————————————————————————————————————————————————————————————
  # Test Core Functionality
  # ——————————————————————————————————————————————————————————————————————————

  cat("\n========================================================================\n")
  cat("  Testing Core Package Loading\n")
  cat("========================================================================\n\n")

  test_packages <- c("DiffBind", "pheatmap", "readr", "dplyr", "ggplot2",
                     "GenomicRanges", "GenomicFeatures")

  all_loaded <- TRUE
  for (pkg in test_packages) {
    cat(sprintf("%-25s ... ", pkg))
    test_result <- tryCatch({
      library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
      cat("✓ OK\n")
      TRUE
    }, error = function(e) {
      cat("✗ FAILED\n")
      FALSE
    })
    if (!test_result) all_loaded <- FALSE
  }

  # ——————————————————————————————————————————————————————————————————————————
  # Summary Report
  # ——————————————————————————————————————————————————————————————————————————

  cat("\n========================================================================\n")
  cat("  Installation Summary\n")
  cat("========================================================================\n\n")

  if (length(installed_pkgs) > 0) {
    cat("✓ Newly installed:", length(installed_pkgs), "packages\n")
  }

  if (length(failed_pkgs) > 0) {
    cat(" Failed:", length(failed_pkgs), "packages:",
        paste(failed_pkgs, collapse = ", "), "\n\n")
  }

  cat("\nR version:", as.character(getRversion()), "\n")
  cat("Bioconductor version:", as.character(BiocManager::version()), "\n")

  # Check key optional packages
  cat("\n--- Optional Features ---\n")
  if (requireNamespace("ChIPseeker", quietly = TRUE)) {
    cat("✓ ChIPseeker: Peak annotation ENABLED\n")
  } else {
    cat(" ChIPseeker: Peak annotation DISABLED\n")
    cat("  Install with: BiocManager::install('ChIPseeker')\n")
  }

  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    cat("✓ rtracklayer: Consensus peaks ENABLED\n")
  } else {
    cat("✗ rtracklayer: Consensus peaks DISABLED\n")
  }

  if (requireNamespace("memes", quietly = TRUE)) {
    cat("✓ memes: R-based motif analysis ENABLED\n")
  } else {
    cat("memes: R-based motif analysis DISABLED (HOMER still works)\n")
  }

  if (requireNamespace("ATACseqQC", quietly = TRUE) &&
      requireNamespace("GenomicAlignments", quietly = TRUE)) {
    cat("✓ ATACseqQC: Quality control metrics ENABLED\n")
  } else {
    cat("ATACseqQC: Quality control metrics DISABLED\n")
    cat("  Install with: BiocManager::install(c('ATACseqQC', 'GenomicAlignments'))\n")
  }

  cat("\n========================================================================\n")

  if (all_loaded) {
    cat("  ✓ SUCCESS! All core packages ready.\n")
    cat("========================================================================\n\n")
    cat("You can now run:\n")
    cat("  Rscript generate_atac_figs_v2.R --base_dir=/path/to/projects\n\n")
  } else {
    cat("   WARNING: Some core packages failed.\n")
    cat("========================================================================\n\n")
    cat("Try updating R and Bioconductor:\n")
    cat("  BiocManager::install(version = '3.18')\n\n")
  }

  cat("Installation complete!\n\n")

  # Cleanup
  rm(cran_packages, bioc_core, bioc_annot, optional, installed_pkgs,
     failed_pkgs, test_packages, all_loaded)



# ——————————————————————————————————————————————————————————————————————————————
# LOAD REQUIRED LIBRARIES
# ——————————————————————————————————————————————————————————————————————————————

suppressPackageStartupMessages({
  library(DiffBind)
  library(pheatmap)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(GenomicRanges)
  library(GenomicFeatures)

  # Optional but recommended packages
  tryCatch({
    library(ChIPseeker)
    CHIPSEEKER_AVAILABLE <- TRUE
  }, error = function(e) {
    warning("ChIPseeker not available. Install with: BiocManager::install('ChIPseeker')")
    CHIPSEEKER_AVAILABLE <<- FALSE
  })

  tryCatch({
    library(rtracklayer)
    RTRACKLAYER_AVAILABLE <- TRUE
  }, error = function(e) {
    warning("rtracklayer not available. Install with: BiocManager::install('rtracklayer')")
    RTRACKLAYER_AVAILABLE <<- FALSE
  })

  tryCatch({
    library(memes)
    MEMES_AVAILABLE <- TRUE
  }, error = function(e) {
    message("memes not available (optional). Install with: BiocManager::install('memes')")
    MEMES_AVAILABLE <<- FALSE
  })
})

# ——————————————————————————————————————————————————————————————————————————————
# CONFIGURATION AND CONSTANTS
# ——————————————————————————————————————————————————————————————————————————————

# Default parameters
DEFAULT_PROMOTER_UPSTREAM <- 3000
DEFAULT_PROMOTER_DOWNSTREAM <- 3000
DEFAULT_TSS_REGION <- c(-3000, 3000)
DEFAULT_MOTIF_WIDTH <- c(6, 20)
DEFAULT_TOP_MOTIFS <- 10

# Figure parameters
FIG_DPI <- 300
FIG_WIDTH <- 7
FIG_HEIGHT <- 7
FIG_WIDTH_WIDE <- 10

# ——————————————————————————————————————————————————————————————————————————————
# UTILITY FUNCTIONS
# ——————————————————————————————————————————————————————————————————————————————

message_step <- function(step_num, msg) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("STEP %d: %s\n", step_num, msg))
  cat(sprintf("========================================\n"))
}

# Check if HOMER is available
check_homer <- function() {
  homer_check <- system("which findMotifsGenome.pl", ignore.stdout = TRUE, ignore.stderr = TRUE)
  return(homer_check == 0)
}

# ——————————————————————————————————————————————————————————————————————————————
# 1. DISCOVER SPECIES FOLDERS
# ——————————————————————————————————————————————————————————————————————————————

find_species_dirs <- function(base_dir) {
  message("→ Searching for species directories in: ", base_dir)

  meta_paths <- list.files(base_dir,
                           pattern    = "^sample_metadata\\.csv$",
                           recursive  = TRUE,
                           full.names = TRUE)

  species_dirs <- dirname(meta_paths)
  message("   Found ", length(species_dirs), " species with metadata files")

  return(species_dirs)
}

# ——————————————————————————————————————————————————————————————————————————————
# 2. READ & SANITIZE METADATA
# ——————————————————————————————————————————————————————————————————————————————

read_metadata <- function(species_dir) {
  peaks_dir <- file.path(species_dir, "peaks_output")

  # Read metadata with required and optional columns
  md <- read_csv(file.path(species_dir, "sample_metadata.csv"),
                 col_types = cols(
                   sampID    = col_character(),
                   condition = col_character(),
                   replicate = col_integer(),
                   TrueID    = col_character(),  # Optional
                   sex       = col_character(),  # Optional: M/F
                   age_category = col_character()  # Optional: Y/O
                 )) %>%
    mutate(
      # Use TrueID for display if present, otherwise sampID
      display_name = ifelse(!is.na(TrueID) & TrueID != "", TrueID, sampID),

      # Validate and standardize sex if present
      sex = if("sex" %in% names(.)) {
        toupper(sex)
      } else {
        NA_character_
      },

      # Validate and standardize age_category if present
      age_category = if("age_category" %in% names(.)) {
        toupper(age_category)
      } else {
        NA_character_
      },

      # File paths based on sampID (folder name)
      peaks_file = normalizePath(
        file.path(peaks_dir, paste0(sampID, "_peaks.narrowPeak")),
        winslash = "/", mustWork = FALSE),
      peaks_xls = normalizePath(
        file.path(peaks_dir, paste0(sampID, "_peaks.xls")),
        winslash = "/", mustWork = FALSE),
      bam_file = normalizePath(
        file.path(peaks_dir, paste0(sampID, "_size_selected_100.bam")),
        winslash = "/", mustWork = FALSE)
    ) %>%
    filter(
      file.exists(peaks_xls),
      !is.na(file.info(peaks_xls)$size),
      file.info(peaks_xls)$size > 0,
      file.exists(bam_file)
    )

  # Report optional column availability
  has_trueid <- any(!is.na(md$TrueID) & md$TrueID != "")
  has_sex <- any(!is.na(md$sex))
  has_age <- any(!is.na(md$age_category))

  if (has_trueid) message("   ✓ TrueID column found - using for display labels")
  if (has_sex) message("   ✓ Sex column found (", paste(unique(na.omit(md$sex)), collapse=", "), ")")
  if (has_age) message("   ✓ Age category found (", paste(unique(na.omit(md$age_category)), collapse=", "), ")")

  return(md)
}

# ——————————————————————————————————————————————————————————————————————————————
# 3. BUILD DIFFBIND OBJECT
# ——————————————————————————————————————————————————————————————————————————————

build_dba <- function(md) {
  message("   Building DiffBind object from ", nrow(md), " samples")

  # Check for optional factors
  has_sex <- "sex" %in% names(md) && any(!is.na(md$sex))
  has_age <- "age_category" %in% names(md) && any(!is.na(md$age_category))

  db <- NULL
  for (i in seq_len(nrow(md))) {
    row <- md[i, ]

    # Base parameters (always included)
    params <- list(
      DBA = db,
      peaks = row$peaks_xls,
      peak.caller = "macs",
      sampID = row$display_name,  # Use TrueID if available
      condition = row$condition,
      replicate = row$replicate,
      bamReads = row$bam_file
    )

    # Add optional factors if present
    if (has_sex && !is.na(row$sex)) {
      params$Treatment <- row$sex  # DiffBind uses "Treatment" for additional factor
    }

    if (has_age && !is.na(row$age_category)) {
      params$Tissue <- row$age_category  # DiffBind uses "Tissue" for another factor
    }

    # Build DBA with dynamic parameters
    db <- do.call(dba.peakset, params)
  }

  if (has_sex) message("   ✓ Sex included as factor")
  if (has_age) message("   ✓ Age category included as factor")

  return(db)
}

# ——————————————————————————————————————————————————————————————————————————————
# 4. ORIGINAL VISUALIZATION FUNCTIONS (Venn, Correlation, Heatmaps)
# ——————————————————————————————————————————————————————————————————————————————

plot_venn_all <- function(db, out_dir, title_prefix) {
  message("   Creating Venn diagram (all samples)")
  png(file.path(out_dir, "venn_all.png"),
      width=FIG_WIDTH, height=FIG_HEIGHT, units="in", res=FIG_DPI)
  dba.plotVenn(db, db$masks$All, main = paste0("All — ", title_prefix))
  dev.off()
}

plot_venn_by_condition <- function(db, md, out_dir, title_prefix) {
  message("   Creating per-condition Venn diagrams")
  samples_df <- db$samples

  for (cond in unique(md$condition)) {
    sel <- which(samples_df$Condition == cond)
    if (length(sel) < 2) {
      message("   • skipping Venn for '", cond, "' (", length(sel), " sample[s])")
      next
    }
    png(file.path(out_dir, paste0("venn_", cond, ".png")),
        width=FIG_WIDTH, height=FIG_HEIGHT, units="in", res=FIG_DPI)
    dba.plotVenn(db, sel, main = paste0(cond, " — ", title_prefix))
    dev.off()
  }
}

plot_venn_by_sex <- function(db, md, out_dir, title_prefix) {
  # Only run if sex column exists
  if (!"sex" %in% names(md) || !any(!is.na(md$sex))) {
    return(NULL)
  }

  message("   Creating per-sex Venn diagrams")
  samples_df <- db$samples

  # DiffBind stores sex in "Treatment" column
  if (!"Treatment" %in% names(samples_df)) {
    return(NULL)
  }

  for (sex_val in unique(na.omit(md$sex))) {
    sel <- which(samples_df$Treatment == sex_val)
    if (length(sel) < 2) {
      message("   • skipping Venn for sex='", sex_val, "' (", length(sel), " sample[s])")
      next
    }
    png(file.path(out_dir, paste0("venn_sex_", sex_val, ".png")),
        width=FIG_WIDTH, height=FIG_HEIGHT, units="in", res=FIG_DPI)
    dba.plotVenn(db, sel, main = paste0("Sex: ", sex_val, " — ", title_prefix))
    dev.off()
  }
}

plot_venn_by_age <- function(db, md, out_dir, title_prefix) {
  # Only run if age_category column exists
  if (!"age_category" %in% names(md) || !any(!is.na(md$age_category))) {
    return(NULL)
  }

  message("   Creating per-age Venn diagrams")
  samples_df <- db$samples

  # DiffBind stores age_category in "Tissue" column
  if (!"Tissue" %in% names(samples_df)) {
    return(NULL)
  }

  for (age_val in unique(na.omit(md$age_category))) {
    sel <- which(samples_df$Tissue == age_val)
    if (length(sel) < 2) {
      message("   • skipping Venn for age='", age_val, "' (", length(sel), " sample[s])")
      next
    }
    png(file.path(out_dir, paste0("venn_age_", age_val, ".png")),
        width=FIG_WIDTH, height=FIG_HEIGHT, units="in", res=FIG_DPI)
    dba.plotVenn(db, sel, main = paste0("Age: ", age_val, " — ", title_prefix))
    dev.off()
  }
}

plot_correlation_and_heatmaps <- function(db, out_dir, title_prefix) {
  message("   Generating correlation plots and heatmaps")

  # Count reads
  dbc <- dba.count(db, bParallel = FALSE, score = DBA_SCORE_READS)
  mat <- dba.peakset(dbc, bRetrieve=TRUE) %>%
    mcols() %>% as.data.frame() %>% as.matrix()
  mode(mat) <- "numeric"

  # Correlation plot
  png(file.path(out_dir, "correlation.png"),
      width=FIG_WIDTH, height=FIG_HEIGHT, units="in", res=FIG_DPI)
  plot(dbc, main = paste0("Correlation — ", title_prefix))
  dev.off()

  # Base heatmap
  png(file.path(out_dir, "heatmap_base.png"),
      width=FIG_WIDTH, height=FIG_HEIGHT, units="in", res=FIG_DPI)
  par(mar=c(5,5,6,2), cex.main=1.2)
  heatmap(mat, main=paste0("Heatmap — ", title_prefix),
          margins=c(5,10), cexRow=0.8, cexCol=0.8)
  dev.off()

  # pheatmap
  png(file.path(out_dir, "heatmap_pheatmap.png"),
      width=FIG_WIDTH+1, height=FIG_HEIGHT+1, units="in", res=FIG_DPI)
  pheatmap(mat,
           main=paste0("Heatmap — ", title_prefix),
           fontsize=10, fontsize_row=8, fontsize_col=8,
           cellwidth=10, cellheight=10,
           treeheight_row=50, treeheight_col=50)
  dev.off()
}

# ——————————————————————————————————————————————————————————————————————————————
# 5. NEW: LOAD GENOME ANNOTATION
# ——————————————————————————————————————————————————————————————————————————————

load_txdb <- function(species_dir) {
  # Try to load TxDb from config or common locations
  config_file <- file.path(species_dir, "config.yaml")

  # Check if GTF/GFF file is specified in config
  if (file.exists(config_file)) {
    config_lines <- readLines(config_file)
    gtf_line <- grep("^GTF_FILE:|^GFF_FILE:", config_lines, value = TRUE)

    if (length(gtf_line) > 0) {
      gtf_path <- sub("^[^:]+:\\s*", "", gtf_line[1])
      gtf_path <- trimws(gtf_path)

      # Handle relative paths
      if (!grepl("^/", gtf_path)) {
        gtf_path <- file.path(species_dir, gtf_path)
      }

      if (file.exists(gtf_path)) {
        message("   Loading annotation from: ", gtf_path)
        txdb <- makeTxDbFromGFF(gtf_path, format = "auto")
        return(txdb)
      }
    }
  }

  # Try to find GTF/GFF in reference directory
  ref_dir <- file.path(species_dir, "reference")
  if (dir.exists(ref_dir)) {
    gtf_files <- list.files(ref_dir, pattern = "\\.(gtf|gff|gff3)(\\.gz)?$",
                           full.names = TRUE, ignore.case = TRUE)
    if (length(gtf_files) > 0) {
      message("   Loading annotation from: ", gtf_files[1])
      txdb <- makeTxDbFromGFF(gtf_files[1], format = "auto")
      return(txdb)
    }
  }

  message("   WARNING: No GTF/GFF annotation file found. Skipping annotation steps.")
  message("   Add GTF_FILE: path/to/annotation.gtf to config.yaml")
  return(NULL)
}

# ——————————————————————————————————————————————————————————————————————————————
# 6. NEW: PEAK FEATURE ANNOTATION
# ——————————————————————————————————————————————————————————————————————————————

annotate_peaks <- function(peaks_file, txdb, out_dir, sample_name) {
  if (is.null(txdb) || !CHIPSEEKER_AVAILABLE) {
    message("   Skipping annotation for ", sample_name, " (no TxDb or ChIPseeker)")
    return(NULL)
  }

  message("   Annotating peaks for: ", sample_name)

  # Read peaks
  peaks <- readPeakFile(peaks_file)

  # Annotate
  peak_anno <- annotatePeak(peaks,
                            tssRegion = DEFAULT_TSS_REGION,
                            TxDb = txdb,
                            annoDb = NULL)  # Can add OrgDb for gene symbols

  # Save annotation table
  anno_df <- as.data.frame(peak_anno)
  write_csv(anno_df, file.path(out_dir, paste0(sample_name, "_peak_annotation.csv")))

  # Generate annotation plots
  png(file.path(out_dir, paste0(sample_name, "_annotation_pie.png")),
      width=FIG_WIDTH, height=FIG_HEIGHT, units="in", res=FIG_DPI)
  print(plotAnnoPie(peak_anno))
  dev.off()

  png(file.path(out_dir, paste0(sample_name, "_annotation_bar.png")),
      width=FIG_WIDTH, height=FIG_HEIGHT, units="in", res=FIG_DPI)
  print(plotAnnoBar(peak_anno))
  dev.off()

  png(file.path(out_dir, paste0(sample_name, "_distance_to_tss.png")),
      width=FIG_WIDTH, height=FIG_HEIGHT, units="in", res=FIG_DPI)
  print(plotDistToTSS(peak_anno, title = paste0("Distance to TSS - ", sample_name)))
  dev.off()

  return(peak_anno)
}

# ——————————————————————————————————————————————————————————————————————————————
# 7. NEW: MOTIF ANALYSIS WITH HOMER
# ——————————————————————————————————————————————————————————————————————————————

run_homer_motif_analysis <- function(peaks_file, genome_fasta, out_dir, sample_name) {
  if (!check_homer()) {
    message("   HOMER not found. Skipping motif analysis for ", sample_name)
    message("   Install HOMER from: http://homer.ucsd.edu/homer/")
    return(FALSE)
  }

  message("   Running HOMER motif analysis for: ", sample_name)

  motif_dir <- file.path(out_dir, paste0(sample_name, "_motifs"))
  dir.create(motif_dir, recursive = TRUE, showWarnings = FALSE)

  # Run HOMER findMotifsGenome.pl
  homer_cmd <- sprintf(
    "findMotifsGenome.pl '%s' '%s' '%s' -size 200 -mask",
    peaks_file,
    genome_fasta,
    motif_dir
  )

  system(homer_cmd)

  # Check if HOMER succeeded
  if (file.exists(file.path(motif_dir, "homerResults.html"))) {
    message("   ✓ HOMER analysis complete. Results in: ", motif_dir)
    return(TRUE)
  } else {
    message("   ✗ HOMER analysis failed")
    return(FALSE)
  }
}

# ——————————————————————————————————————————————————————————————————————————————
# 8. NEW: COMPREHENSIVE PEAK FEATURE TABLE
# ——————————————————————————————————————————————————————————————————————————————

create_peak_features_table <- function(md, species_dir, out_dir) {
  message("   Creating comprehensive peak features table")

  all_features <- list()

  for (i in seq_len(nrow(md))) {
    row <- md[i, ]
    sample_name <- row$display_name  # Use TrueID if available
    peaks_file <- row$peaks_file

    if (!file.exists(peaks_file)) {
      message("   WARNING: Peak file not found for ", sample_name)
      next
    }

    # Read narrowPeak file
    peaks <- read.table(peaks_file,
                       col.names = c("chr", "start", "end", "name", "score",
                                   "strand", "signalValue", "pValue", "qValue", "peak"))

    # Calculate features (base columns)
    features <- data.frame(
      sample = sample_name,
      condition = row$condition,
      replicate = row$replicate,
      total_peaks = nrow(peaks),
      mean_peak_width = mean(peaks$end - peaks$start),
      median_peak_width = median(peaks$end - peaks$start),
      mean_signal = mean(peaks$signalValue),
      median_signal = median(peaks$signalValue),
      mean_qvalue = mean(peaks$qValue),
      peaks_high_confidence = sum(peaks$qValue > 10),  # q > 10
      peaks_strong_signal = sum(peaks$signalValue > quantile(peaks$signalValue, 0.75)),
      stringsAsFactors = FALSE
    )

    # Add optional metadata columns
    if ("sex" %in% names(row) && !is.na(row$sex)) {
      features$sex <- row$sex
    }
    if ("age_category" %in% names(row) && !is.na(row$age_category)) {
      features$age_category <- row$age_category
    }

    all_features[[sample_name]] <- features
  }

  # Combine all features
  features_table <- bind_rows(all_features)

  # Save table
  write_csv(features_table, file.path(out_dir, "peak_features_summary.csv"))

  # Create visualization
  if (nrow(features_table) > 0) {
    png(file.path(out_dir, "peak_counts_by_sample.png"),
        width=FIG_WIDTH_WIDE, height=FIG_HEIGHT, units="in", res=FIG_DPI)

    p <- ggplot(features_table, aes(x = sample, y = total_peaks, fill = condition)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Total Peaks per Sample",
           x = "Sample", y = "Number of Peaks")
    print(p)
    dev.off()

    # Peak width distribution
    png(file.path(out_dir, "peak_width_distribution.png"),
        width=FIG_WIDTH_WIDE, height=FIG_HEIGHT, units="in", res=FIG_DPI)

    p <- ggplot(features_table, aes(x = sample, y = mean_peak_width, fill = condition)) +
      geom_bar(stat = "identity") +
      geom_errorbar(aes(ymin = median_peak_width, ymax = mean_peak_width), width = 0.2) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Peak Width Distribution",
           x = "Sample", y = "Mean Peak Width (bp)")
    print(p)
    dev.off()
  }

  return(features_table)
}

# ——————————————————————————————————————————————————————————————————————————————
# 9. NEW: GENOMIC DISTRIBUTION ANALYSIS
# ——————————————————————————————————————————————————————————————————————————————

analyze_genomic_distribution <- function(peak_annotations, out_dir, species_name) {
  if (length(peak_annotations) == 0 || !CHIPSEEKER_AVAILABLE) {
    return(NULL)
  }

  message("   Analyzing genomic distribution across all samples")

  # Combine all annotations
  png(file.path(out_dir, "genomic_distribution_comparison.png"),
      width=FIG_WIDTH_WIDE, height=FIG_HEIGHT, units="in", res=FIG_DPI)
  print(plotAnnoBar(peak_annotations))
  dev.off()

  # Distance to TSS comparison
  png(file.path(out_dir, "distance_to_tss_comparison.png"),
      width=FIG_WIDTH_WIDE, height=FIG_HEIGHT, units="in", res=FIG_DPI)
  print(plotDistToTSS(peak_annotations))
  dev.off()

  # Create summary table
  distribution_summary <- lapply(names(peak_annotations), function(name) {
    anno <- as.data.frame(peak_annotations[[name]])

    data.frame(
      sample = name,
      promoter = sum(grepl("Promoter", anno$annotation)) / nrow(anno) * 100,
      exon = sum(grepl("Exon", anno$annotation)) / nrow(anno) * 100,
      intron = sum(grepl("Intron", anno$annotation)) / nrow(anno) * 100,
      intergenic = sum(grepl("Intergenic", anno$annotation)) / nrow(anno) * 100,
      mean_dist_to_tss = mean(abs(anno$distanceToTSS), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  distribution_df <- bind_rows(distribution_summary)
  write_csv(distribution_df, file.path(out_dir, "genomic_distribution_summary.csv"))

  return(distribution_df)
}

# ——————————————————————————————————————————————————————————————————————————————
# 10. NEW: CONSENSUS PEAK SET
# ——————————————————————————————————————————————————————————————————————————————

create_consensus_peaks <- function(md, out_dir, condition = NULL) {
  message("   Creating consensus peak set")

  # Filter by condition if specified
  if (!is.null(condition)) {
    md <- md %>% filter(condition == !!condition)
    suffix <- paste0("_", condition)
  } else {
    suffix <- "_all"
  }

  if (nrow(md) == 0) return(NULL)

  # Read all peaks
  all_peaks <- lapply(md$peaks_file, function(f) {
    if (file.exists(f)) {
      rtracklayer::import(f, format = "bed")
    } else {
      NULL
    }
  })

  all_peaks <- all_peaks[!sapply(all_peaks, is.null)]

  if (length(all_peaks) == 0) return(NULL)

  # Combine and reduce to consensus
  combined <- do.call(c, all_peaks)
  consensus <- reduce(combined)

  # Save consensus peaks
  consensus_file <- file.path(out_dir, paste0("consensus_peaks", suffix, ".bed"))
  rtracklayer::export(consensus, consensus_file, format = "bed")

  message("   ✓ Consensus peaks saved: ", consensus_file, " (", length(consensus), " peaks)")

  return(consensus)
}

# ——————————————————————————————————————————————————————————————————————————————
# 11. NEW: ATAC-seq QUALITY CONTROL METRICS
# ——————————————————————————————————————————————————————————————————————————————

calculate_complexity_metrics <- function(bam_file) {
  message("   Calculating complexity metrics for: ", basename(bam_file))

  # Read alignments
  reads <- GenomicAlignments::readGAlignments(bam_file)
  total_reads <- length(reads)

  # Create unique position identifiers
  positions <- paste(seqnames(reads), start(reads), end(reads), strand(reads), sep="_")
  unique_positions <- length(unique(positions))

  # Calculate NRF (Non-Redundant Fraction)
  NRF <- unique_positions / total_reads

  # Calculate duplication statistics
  pos_counts <- table(positions)
  dup_freq_table <- table(pos_counts)

  # M1 = positions with exactly 1 read
  M1 <- as.numeric(dup_freq_table["1"])
  if(is.na(M1)) M1 <- 0

  # M2 = positions with exactly 2 reads
  M2 <- as.numeric(dup_freq_table["2"])
  if(is.na(M2)) M2 <- 0

  # Calculate PBC metrics
  PBC1 <- M1 / unique_positions
  PBC2 <- ifelse(M2 > 0, M1 / M2, NA)

  # Duplication rate
  dup_rate <- 1 - NRF

  # Quality assessment
  nrf_quality <- ifelse(NRF > 0.9, "Excellent",
                        ifelse(NRF > 0.8, "Good",
                               ifelse(NRF > 0.5, "Acceptable", "Concerning")))

  pbc1_quality <- ifelse(PBC1 > 0.9, "Excellent",
                         ifelse(PBC1 > 0.7, "Good",
                                ifelse(PBC1 > 0.5, "Acceptable", "Concerning")))

  # Overall quality
  if(NRF > 0.8 && PBC1 > 0.7) {
    overall_quality <- "EXCELLENT"
  } else if(NRF > 0.5 && PBC1 > 0.5) {
    overall_quality <- "ACCEPTABLE"
  } else {
    overall_quality <- "CONCERNING"
  }

  return(list(
    total_reads = total_reads,
    unique_positions = unique_positions,
    NRF = NRF,
    PBC1 = PBC1,
    PBC2 = PBC2,
    M1 = M1,
    M2 = M2,
    dup_rate = dup_rate,
    nrf_quality = nrf_quality,
    pbc1_quality = pbc1_quality,
    overall_quality = overall_quality
  ))
}

analyze_fragment_sizes <- function(bam_file) {
  message("   Analyzing fragment sizes for: ", basename(bam_file))

  tryCatch({
    # Read as pairs
    pairs <- GenomicAlignments::readGAlignmentPairs(bam_file)
    n_pairs <- length(pairs)

    # Calculate insert sizes
    first_read <- GenomicAlignments::first(pairs)
    last_read <- GenomicAlignments::last(pairs)
    insert_sizes <- abs(end(last_read) - start(first_read)) + 1

    # Filter for typical ATAC-seq range
    insert_sizes_filtered <- insert_sizes[insert_sizes <= 1000 & insert_sizes > 0]

    # Calculate statistics
    median_frag <- median(insert_sizes_filtered)
    mean_frag <- mean(insert_sizes_filtered)

    # Count fragments in different size ranges
    nfr_count <- sum(insert_sizes_filtered < 100)
    mono_count <- sum(insert_sizes_filtered >= 100 & insert_sizes_filtered < 300)
    di_count <- sum(insert_sizes_filtered >= 300 & insert_sizes_filtered < 500)

    nfr_pct <- (nfr_count / length(insert_sizes_filtered)) * 100
    mono_pct <- (mono_count / length(insert_sizes_filtered)) * 100
    di_pct <- (di_count / length(insert_sizes_filtered)) * 100

    # Create histogram
    frag_hist <- hist(insert_sizes_filtered,
                      breaks=seq(0, 1000, by=10),
                      plot=FALSE)

    return(list(
      n_pairs = n_pairs,
      median = median_frag,
      mean = mean_frag,
      sizes = insert_sizes_filtered,
      hist = frag_hist,
      nfr_count = nfr_count,
      mono_count = mono_count,
      di_count = di_count,
      nfr_pct = nfr_pct,
      mono_pct = mono_pct,
      di_pct = di_pct
    ))
  }, error = function(e) {
    message("   ERROR: ", e$message)
    return(NULL)
  })
}

create_qc_plots <- function(complexity_metrics, fragment_data, sample_name, out_dir) {
  message("   Creating QC plots for: ", sample_name)

  # Extract metrics
  NRF <- complexity_metrics$NRF
  PBC1 <- complexity_metrics$PBC1

  # Plot 1: Library Complexity Metrics
  png(file.path(out_dir, paste0("QC_metrics_", sample_name, ".png")),
      width=8, height=6, units="in", res=300)

  par(mfrow=c(1,2), mar=c(5,5,4,2))

  # NRF barplot
  barplot(NRF, ylim=c(0,1),
          main=paste("Non-Redundant Fraction\n", sample_name),
          ylab="NRF", col="white", border="black", names.arg="")
  abline(h=c(0.9, 0.8, 0.5), lty=c(2,3,4), col="black", lwd=1.5)
  text(0.7, 0.92, "Excellent (>0.9)", cex=0.8)
  text(0.7, 0.82, "Good (>0.8)", cex=0.8)
  text(0.7, NRF - 0.05, sprintf("%.4f", NRF), cex=1.2, font=2)

  # PBC1 barplot
  barplot(PBC1, ylim=c(0,1),
          main=paste("PCR Bottleneck Coefficient\n", sample_name),
          ylab="PBC1", col="white", border="black", names.arg="")
  abline(h=c(0.9, 0.7, 0.5), lty=c(2,3,4), col="black", lwd=1.5)
  text(0.7, 0.92, "Excellent (>0.9)", cex=0.8)
  text(0.7, 0.72, "Good (>0.7)", cex=0.8)
  text(0.7, PBC1 - 0.05, sprintf("%.4f", PBC1), cex=1.2, font=2)

  dev.off()

  # Plot 2: Fragment Size Distribution
  if(!is.null(fragment_data)) {
    png(file.path(out_dir, paste0("fragsize_", sample_name, ".png")),
        width=10, height=6, units="in", res=300)

    par(mar=c(5,5,4,2))
    plot(fragment_data$hist$mids, fragment_data$hist$counts,
         type="l", lwd=2, col="black",
         main=paste("Fragment Size Distribution -", sample_name),
         xlab="Fragment Length (bp)",
         ylab="Number of Fragments",
         xlim=c(0, 1000))

    # Add reference lines
    abline(v=c(100, 200, 400), lty=2, col="gray30", lwd=1.5)

    # Add labels
    text(50, max(fragment_data$hist$counts)*0.9, "NFR\n(<100bp)", cex=0.9)
    text(150, max(fragment_data$hist$counts)*0.7, "Mono", cex=0.9)
    text(350, max(fragment_data$hist$counts)*0.5, "Di", cex=0.9)

    dev.off()

    # Plot 3: Combined QC Overview
    png(file.path(out_dir, paste0("combined_QC_", sample_name, ".png")),
        width=12, height=8, units="in", res=300)

    layout(matrix(c(1,1,2,
                    1,1,3), nrow=2, byrow=TRUE))

    # Fragment size (large, left)
    par(mar=c(5,5,4,2))
    plot(fragment_data$hist$mids, fragment_data$hist$counts,
         type="l", lwd=2, col="black",
         main=paste("ATAC-seq QC Summary -", sample_name),
         xlab="Fragment Length (bp)",
         ylab="Number of Fragments",
         xlim=c(0, 1000))
    abline(v=c(100, 200, 400), lty=2, col="gray30")

    # NRF (small, top right)
    par(mar=c(4,4,3,1))
    barplot(NRF, ylim=c(0,1),
            main="NRF", ylab="", col="white", border="black")
    abline(h=c(0.8, 0.9), lty=2, col="gray30")
    text(0.7, NRF-0.05, sprintf("%.3f", NRF), font=2)

    # PBC1 (small, bottom right)
    par(mar=c(4,4,3,1))
    barplot(PBC1, ylim=c(0,1),
            main="PBC1", ylab="", col="white", border="black")
    abline(h=c(0.7, 0.9), lty=2, col="gray30")
    text(0.7, PBC1-0.05, sprintf("%.3f", PBC1), font=2)

    dev.off()
  }

  return(list(
    metrics_plot = file.path(out_dir, paste0("QC_metrics_", sample_name, ".png")),
    fragsize_plot = file.path(out_dir, paste0("fragsize_", sample_name, ".png")),
    combined_plot = file.path(out_dir, paste0("combined_QC_", sample_name, ".png"))
  ))
}

run_atacseq_qc <- function(md, species_dir, out_dir) {
  message("   Running ATAC-seq Quality Control")

  # Find BAM files
  bam_files <- md$bam_file[file.exists(md$bam_file)]

  if(length(bam_files) == 0) {
    message("   WARNING: No BAM files found")
    return(NULL)
  }

  message("   Found ", length(bam_files), " BAM files")

  # Initialize results
  all_results <- list()

  # Process each BAM file
  for(i in seq_along(bam_files)) {
    bam_file <- bam_files[i]
    sample_name <- md$display_name[md$bam_file == bam_file][1]

    message("\n   Processing: ", sample_name)

    # Calculate metrics
    complexity <- calculate_complexity_metrics(bam_file)
    fragment_data <- analyze_fragment_sizes(bam_file)

    # Create plots
    plots <- create_qc_plots(complexity, fragment_data, sample_name, out_dir)

    # Store results
    all_results[[sample_name]] <- list(
      complexity = complexity,
      fragments = fragment_data,
      plots = plots
    )
  }

  # Create summary table
  summary_df <- data.frame(
    Sample = names(all_results),
    Total_Reads = sapply(all_results, function(x) x$complexity$total_reads),
    Unique_Pos = sapply(all_results, function(x) x$complexity$unique_positions),
    NRF = sapply(all_results, function(x) round(x$complexity$NRF, 4)),
    PBC1 = sapply(all_results, function(x) round(x$complexity$PBC1, 4)),
    PBC2 = sapply(all_results, function(x) round(x$complexity$PBC2, 4)),
    Dup_Rate_Pct = sapply(all_results, function(x) round(x$complexity$dup_rate * 100, 2)),
    NRF_Quality = sapply(all_results, function(x) x$complexity$nrf_quality),
    PBC1_Quality = sapply(all_results, function(x) x$complexity$pbc1_quality),
    Overall = sapply(all_results, function(x) x$complexity$overall_quality),
    stringsAsFactors = FALSE
  )

  # Add fragment size info if available
  if(!is.null(all_results[[1]]$fragments)) {
    summary_df$Median_FragSize <- sapply(all_results, function(x) {
      ifelse(!is.null(x$fragments), x$fragments$median, NA)
    })
    summary_df$NFR_Pct <- sapply(all_results, function(x) {
      ifelse(!is.null(x$fragments), round(x$fragments$nfr_pct, 1), NA)
    })
    summary_df$Mono_Pct <- sapply(all_results, function(x) {
      ifelse(!is.null(x$fragments), round(x$fragments$mono_pct, 1), NA)
    })
  }

  message("\n   QC Summary:")
  print(summary_df)

  return(summary_df)
}

# ——————————————————————————————————————————————————————————————————————————————
# 12. MASTER WRAPPER FUNCTION
# ——————————————————————————————————————————————————————————————————————————————

generate_comprehensive_analysis <- function(base_dir,
                                           run_motifs = TRUE,
                                           run_annotation = TRUE) {

  cat("\n")
  cat("========================================================================\n")
  cat("  IISAGE ATAC-seq Comprehensive Analysis v2.0\n")
  cat("========================================================================\n")
  cat("Base directory: ", base_dir, "\n")
  cat("Run motif analysis: ", run_motifs, "\n")
  cat("Run peak annotation: ", run_annotation, "\n")
  cat("\n")

  species_dirs <- find_species_dirs(base_dir)

  if (length(species_dirs) == 0) {
    stop("No species directories found with sample_metadata.csv")
  }

  for (species_dir in species_dirs) {
    species_name <- basename(species_dir)

    message_step(0, paste("Processing Species:", species_name))
    message("→ Species directory: ", species_dir)

    # Read metadata
    message_step(1, "Reading Metadata")
    md <- read_metadata(species_dir)

    if (nrow(md) == 0) {
      warning("   • no valid samples, skipping ", species_name, "\n")
      next
    }
    message("   ✓ Found ", nrow(md), " valid samples")

    # Create output directory
    out_dir <- file.path(species_dir, "figures")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Build DiffBind object
    message_step(2, "Building DiffBind Object")
    db <- build_dba(md)
    message("   ✓ DiffBind object created")

    # Original visualizations
    message_step(3, "Generating Venn Diagrams")
    plot_venn_all(db, out_dir, species_name)
    plot_venn_by_condition(db, md, out_dir, species_name)

    # Stratified Venn diagrams (if sex/age present)
    plot_venn_by_sex(db, md, out_dir, species_name)
    plot_venn_by_age(db, md, out_dir, species_name)

    message_step(4, "Generating Correlation Plots and Heatmaps")
    plot_correlation_and_heatmaps(db, out_dir, species_name)

    # Load genome annotation
    message_step(5, "Loading Genome Annotation")
    txdb <- NULL
    if (run_annotation) {
      txdb <- load_txdb(species_dir)
    }

    # Peak annotation
    peak_annotations <- list()
    if (run_annotation && !is.null(txdb) && CHIPSEEKER_AVAILABLE) {
      message_step(6, "Annotating Peaks")

      for (i in seq_len(nrow(md))) {
        row <- md[i, ]
        if (file.exists(row$peaks_file)) {
          anno <- annotate_peaks(row$peaks_file, txdb, out_dir, row$display_name)
          if (!is.null(anno)) {
            peak_annotations[[row$display_name]] <- anno
          }
        }
      }

      if (length(peak_annotations) > 0) {
        message_step(7, "Analyzing Genomic Distribution")
        analyze_genomic_distribution(peak_annotations, out_dir, species_name)
      }
    } else {
      message("   Skipping peak annotation (TxDb not available or ChIPseeker not installed)")
    }

    # Create peak features table
    message_step(8, "Creating Peak Features Table")
    features_table <- create_peak_features_table(md, species_dir, out_dir)
    message("   ✓ Feature table saved")

    # Consensus peaks
    if (RTRACKLAYER_AVAILABLE) {
      message_step(9, "Creating Consensus Peak Sets")
      consensus_all <- create_consensus_peaks(md, out_dir)

      # Per-condition consensus
      for (cond in unique(md$condition)) {
        consensus_cond <- create_consensus_peaks(md, out_dir, condition = cond)
      }
    } else {
      message("   Skipping consensus peaks (rtracklayer not available)")
    }

    # Motif analysis
    if (run_motifs) {
      message_step(10, "Running Motif Analysis")

      # Get genome FASTA
      config_file <- file.path(species_dir, "config.yaml")
      genome_fasta <- NULL

      if (file.exists(config_file)) {
        config_lines <- readLines(config_file)
        fasta_line <- grep("^REF_GENOME_FASTA:", config_lines, value = TRUE)

        if (length(fasta_line) > 0) {
          genome_fasta <- sub("^[^:]+:\\s*", "", fasta_line[1])
          genome_fasta <- trimws(genome_fasta)

          if (!grepl("^/", genome_fasta)) {
            genome_fasta <- file.path(species_dir, genome_fasta)
          }
        }
      }

      if (!is.null(genome_fasta) && file.exists(genome_fasta)) {
        for (i in seq_len(nrow(md))) {
          row <- md[i, ]
          if (file.exists(row$peaks_file)) {
            run_homer_motif_analysis(row$peaks_file, genome_fasta, out_dir, row$display_name)
          }
        }
      } else {
        message("   WARNING: Genome FASTA not found. Skipping motif analysis.")
        message("   Add REF_GENOME_FASTA: path/to/genome.fasta to config.yaml")
      }
    }

    # Quality Control Metrics
    if (requireNamespace("ATACseqQC", quietly = TRUE) &&
        requireNamespace("GenomicAlignments", quietly = TRUE) &&
        requireNamespace("Rsamtools", quietly = TRUE)) {
      message_step(11, "Running ATAC-seq Quality Control")

      qc_results <- run_atacseq_qc(md, species_dir, out_dir)

      if (!is.null(qc_results)) {
        write_csv(qc_results, file.path(out_dir, "QC_summary.csv"))
        message("   ✓ QC summary saved: QC_summary.csv")
      }
    } else {
      message("   Skipping QC (ATACseqQC/GenomicAlignments not available)")
      message("   Install with: BiocManager::install(c('ATACseqQC', 'GenomicAlignments'))")
    }

    # Summary
    message_step(12, "Analysis Complete for", species_name)
    message("   Output directory: ", out_dir)
    message("   Total samples processed: ", nrow(md))
    message("   Figures generated: Venn diagrams, correlations, heatmaps")
    if (length(peak_annotations) > 0) {
      message("   Peak annotations: ", length(peak_annotations), " samples")
    }
    if (run_motifs && check_homer()) {
      message("   Motif analysis: HOMER results in *_motifs/ subdirectories")
    }
    if (exists("qc_results") && !is.null(qc_results)) {
      message("   QC metrics: Library complexity, fragment sizes, QC plots")
    }

    cat("\n")
  }

  cat("========================================================================\n")
  cat("  All species processed successfully!\n")
  cat("========================================================================\n\n")
}

# ——————————————————————————————————————————————————————————————————————————————
# COMMAND LINE INTERFACE
# ——————————————————————————————————————————————————————————————————————————————

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0) {
    cat("Usage: Rscript generate_atac_figs_v2.R --base_dir=/path/to/projects [--no-motifs] [--no-annotation]\n")
    quit(status = 1)
  }

  base_dir <- NULL
  run_motifs <- TRUE
  run_annotation <- TRUE

  for (arg in args) {
    if (grepl("^--base_dir=", arg)) {
      base_dir <- sub("^--base_dir=", "", arg)
    } else if (arg == "--no-motifs") {
      run_motifs <- FALSE
    } else if (arg == "--no-annotation") {
      run_annotation <- FALSE
    }
  }

  if (is.null(base_dir)) {
    stop("--base_dir is required")
  }

  generate_comprehensive_analysis(base_dir,
                                  run_motifs = run_motifs,
                                  run_annotation = run_annotation)
}

# ——————————————————————————————————————————————————————————————————————————————
# Example usage from R:
# source("generate_atac_figs_v2.R")
# generate_comprehensive_analysis("/path/to/atac_seq_projects")
# ——————————————————————————————————————————————————————————————————————————————

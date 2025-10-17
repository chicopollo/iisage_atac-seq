library(DiffBind)
library(pheatmap)
library(readr)
library(dplyr)

# ——————————————————————————————————————————————————————————————
# 1. Discover species folders by metadata CSV
# ——————————————————————————————————————————————————————————————
find_species_dirs <- function(base_dir) {
  meta_paths <- list.files(base_dir,
                           pattern    = "^sample_metadata\\.csv$",
                           recursive  = TRUE,
                           full.names = TRUE)
  dirname(meta_paths)
}
species_dir <-find_species_dirs("atac_sec/")
# ——————————————————————————————————————————————————————————————
# 2. Read & sanitize metadata, build file paths
# ——————————————————————————————————————————————————————————————
read_metadata <- function(species_dir) {
  peaks_dir <- file.path(species_dir, "peaks_output")
  md <- read_csv(file.path(species_dir, "sample_metadata.csv"),
                 col_types = cols(
                   sampID    = col_character(),
                   condition = col_character(),
                   replicate = col_integer()
                 )) %>%
    mutate(
      peaks_file = normalizePath(
        file.path(peaks_dir, paste0(sampID, "_peaks.xls")),
        winslash = "/", mustWork = FALSE),
      bam_file = normalizePath(
        file.path(peaks_dir, paste0(sampID, "_size_selected_100.bam")),
        winslash = "/", mustWork = FALSE)
    ) %>%
    filter(
      file.exists(peaks_file),
      !is.na(file.info(peaks_file)$size),
      file.info(peaks_file)$size > 0,
      file.exists(bam_file)
    )
  md
}
md <- read_metadata(species_dir)

# ——————————————————————————————————————————————————————————————
# 3. Build a DiffBind object from metadata
# ——————————————————————————————————————————————————————————————
build_dba <- function(md) {
  db <- NULL
  for (i in seq_len(nrow(md))) {
    row <- md[i, ]
    db <- dba.peakset(db,
                      peaks       = row$peaks_file,
                      peak.caller = "macs",
                      sampID      = row$sampID,
                      condition   = row$condition,
                      replicate   = row$replicate,
                      bamReads    = row$bam_file)
  }
  db
}
build_dba(md)

# ——————————————————————————————————————————————————————————————
# 4a. Plot “All” Venn
# ——————————————————————————————————————————————————————————————
plot_venn_all <- function(db, out_dir, title_prefix) {
  png(file.path(out_dir, "venn_all.png"), width=7, height=7, units="in", res=300)
  dba.plotVenn(db, db$masks$All, main = paste0("All — ", title_prefix))
  dev.off()
}

# ——————————————————————————————————————————————————————————————
# 4b. Plot per-condition Venn
#     (by selecting sample indices)
# ——————————————————————————————————————————————————————————————
plot_venn_by_condition <- function(db, md, out_dir, title_prefix) {
  samples_df <- db$samples
  for (cond in unique(md$condition)) {
    sel <- which(samples_df$Condition == cond)
    if (length(sel) < 2) {
      message("   • skipping Venn for ‘", cond, "’ (", length(sel), " sample[s])")
      next
    }
    png(file.path(out_dir, paste0("venn_", cond, ".png")),
        width=7, height=7, units="in", res=300)
    dba.plotVenn(db, sel, main = paste0(cond, " — ", title_prefix))
    dev.off()
  }
}

# ——————————————————————————————————————————————————————————————
# 5. Correlation plot + heatmaps
# ——————————————————————————————————————————————————————————————
plot_correlation_and_heatmaps <- function(db, out_dir, title_prefix) {
  # counts
  dbc <- dba.count(db, bParallel = FALSE, score = DBA_SCORE_READS)
  mat <- dba.peakset(dbc, bRetrieve=TRUE) %>%
    mcols() %>% as.data.frame() %>% as.matrix()
  mode(mat) <- "numeric"

  # 5a. correlation
  png(file.path(out_dir, "correlation.png"), width=7, height=7, units="in", res=300)
  plot(dbc, main = paste0("Correlation — ", title_prefix))
  dev.off()

  # 5b. base heatmap
  png(file.path(out_dir, "heatmap_base.png"), width=7, height=7, units="in", res=300)
  par(mar=c(5,5,6,2), cex.main=1.2)
  heatmap(mat, main=paste0("Heatmap — ", title_prefix),
          margins=c(5,10), cexRow=0.8, cexCol=0.8)
  dev.off()

  # 5c. pheatmap
  png(file.path(out_dir, "heatmap_pheatmap.png"), width=8, height=8, units="in", res=300)
  pheatmap(mat,
           main=paste0("Heatmap — ", title_prefix),
           fontsize=10, fontsize_row=8, fontsize_col=8,
           cellwidth=10, cellheight=10,
           treeheight_row=50, treeheight_col=50)
  dev.off()
}

# ——————————————————————————————————————————————————————————————
# 6. Master wrapper
# ——————————————————————————————————————————————————————————————
generate_figures_modular <- function(base_dir) {
  for (species_dir in find_species_dirs(base_dir)) {
    message("→ Processing: ", species_dir)
    md      <- read_metadata(species_dir)
    if (nrow(md) == 0) {
      warning("   • no valid samples, skipping\n")
      next
    }
    db      <- build_dba(md)
    out_dir <- file.path(species_dir, "figures")
    dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

    plot_venn_all(db, out_dir, basename(species_dir))
    plot_venn_by_condition(db, md, out_dir, basename(species_dir))
    plot_correlation_and_heatmaps(db, out_dir, basename(species_dir))
  }
}

# ——————————————————————————————————————————————————————————————
# Example:
# generate_figures_modular("atac_sec/")
# —————————————————————————————————————————————

# R Analysis Integration Guide

## Overview

The IISAGE ATAC-seq pipeline v2.0 now includes **automatic R-based downstream analysis** as Stage 4 of the pipeline. This integration provides "one-click" access to comprehensive peak analysis, visualization, and motif discovery.

## Quick Start

### Prerequisites

1. **Create sample_metadata.csv** in your species directory:

```csv
sampID,condition,replicate
sample1,control,1
sample2,control,2
sample3,treatment,1
sample4,treatment,2
```

**Important:**
- `sampID` must match your peak file names (e.g., `sample1_peaks.narrowPeak`)
- Include all samples you want analyzed
- Use consistent condition names

2. **Optional:** Add GTF/GFF annotation file to config.yaml:

```yaml
# In your species/config.yaml
GTF_FILE: reference/annotation.gtf
```

### Basic Usage

```bash
# Standard pipeline with R analysis (default)
sbatch --export=SPECIES_NAME=my_species multi_species_atac_pipeline_v2.sh
```

That's it! R analysis runs automatically after peak calling.

---

## Configuration Options

### Enable/Disable R Analysis

```bash
# Disable R analysis entirely
sbatch --export=SPECIES_NAME=my_species,RUN_R_ANALYSIS=false multi_species_atac_pipeline_v2.sh

# Enable R analysis (default)
sbatch --export=SPECIES_NAME=my_species,RUN_R_ANALYSIS=true multi_species_atac_pipeline_v2.sh
```

### Skip Specific Analyses

```bash
# Skip motif analysis (faster, useful if HOMER not installed)
sbatch --export=SPECIES_NAME=my_species,R_SKIP_MOTIFS=true multi_species_atac_pipeline_v2.sh

# Skip peak annotation (useful if no GTF available)
sbatch --export=SPECIES_NAME=my_species,R_SKIP_ANNOTATION=true multi_species_atac_pipeline_v2.sh

# Skip both
sbatch --export=SPECIES_NAME=my_species,R_SKIP_MOTIFS=true,R_SKIP_ANNOTATION=true multi_species_atac_pipeline_v2.sh
```

### Require R Analysis

By default, pipeline continues even if R analysis fails. To make R analysis mandatory:

```bash
sbatch --export=SPECIES_NAME=my_species,R_ANALYSIS_REQUIRED=true multi_species_atac_pipeline_v2.sh
```

### Custom R Module

If your HPCC uses a different R module name:

```bash
sbatch --export=SPECIES_NAME=my_species,R_MODULE=R/4.4.0 multi_species_atac_pipeline_v2.sh
```

---

## Output Files

R analysis creates a `figures/` directory in your species folder:

```
species_name/
└── figures/
    ├── venn_all.png                           # Venn diagram (all samples)
    ├── venn_control.png                       # Per-condition Venn
    ├── venn_treatment.png
    ├── correlation.png                        # Sample correlation
    ├── heatmap_pheatmap.png                   # Peak intensity heatmap
    ├── genomic_distribution_comparison.png    # Feature distribution
    ├── peak_features_summary.csv              # Peak statistics table
    ├── genomic_distribution_summary.csv       # Genomic context table
    ├── consensus_peaks_all.bed                # Consensus peak set
    ├── sample1_peak_annotation.csv            # Per-sample annotations
    ├── sample1_annotation_pie.png
    ├── sample1_distance_to_tss.png
    └── sample1_motifs/                        # HOMER motif results
        ├── homerResults.html
        └── knownResults.html
```

---

## Usage Scenarios

### Scenario 1: Full Analysis (Default)

```bash
# Run complete pipeline including R analysis
sbatch --export=SPECIES_NAME=chrysemys_picta multi_species_atac_pipeline_v2.sh
```

**What runs:**
- Trimming & alignment
- Peak calling
- DeepTools visualization
- **R analysis** (Venn, annotation, motifs, features)

### Scenario 2: Array Job with R Analysis

```bash
# Parallel sample processing + R analysis
sbatch --array=1-12 --export=SPECIES_NAME=chrysemys_picta multi_species_atac_pipeline_v2.sh
```

**What happens:**
- Array tasks (1-12): Process samples in parallel
- Main job (no array ID): Runs R analysis after all samples complete

### Scenario 3: Quick Analysis (No Motifs)

```bash
# Skip time-consuming motif discovery
sbatch --export=SPECIES_NAME=chrysemys_picta,R_SKIP_MOTIFS=true multi_species_atac_pipeline_v2.sh
```

**Use when:**
- Testing pipeline
- HOMER not installed
- Quick results needed

### Scenario 4: Rerun R Analysis Only

```bash
# Only run R analysis (skip alignment/peaks)
sbatch --export=SPECIES_NAME=chrysemys_picta,RUN_TRIMMING=false,RUN_ALIGNMENT=false,RUN_PEAK_CALLING=false,RUN_DEEPTOOLS=false multi_species_atac_pipeline_v2.sh
```

**Use when:**
- Peak files already exist
- Want to regenerate figures with different settings
- Changed sample_metadata.csv

### Scenario 5: Without R Analysis

```bash
# Traditional pipeline only
sbatch --export=SPECIES_NAME=chrysemys_picta,RUN_R_ANALYSIS=false multi_species_atac_pipeline_v2.sh
```

**Use when:**
- Running R analysis separately later
- No R environment available
- Only need peak calls

---

## Troubleshooting

### Error: "sample_metadata.csv not found"

**Solution:** Create the file in your species directory:

```bash
cd /path/to/species_dir
cat > sample_metadata.csv << EOF
sampID,condition,replicate
sample1,control,1
sample2,control,2
EOF
```

### Error: "Failed to load R module"

**Cause:** R module name doesn't match your HPCC

**Solution:** Find available R modules:

```bash
module spider R
```

Then specify the correct module:

```bash
sbatch --export=SPECIES_NAME=my_species,R_MODULE=R/4.4.1 multi_species_atac_pipeline_v2.sh
```

### Error: "R script not found"

**Cause:** `generate_atac_figs_v2.R` not in scripts directory

**Solution:** Ensure script is in same directory as pipeline:

```bash
ls /path/to/scripts/generate_atac_figs_v2.R
```

Or specify custom path:

```bash
export R_SCRIPT_PATH=/custom/path/generate_atac_figs_v2.R
sbatch --export=SPECIES_NAME=my_species multi_species_atac_pipeline_v2.sh
```

### Warning: "Only N samples (minimum 2 recommended)"

**Cause:** Too few samples for meaningful comparative analysis

**Impact:** R analysis will run but statistical power is limited

**Solution:** Include more biological replicates if possible

### R Analysis Failed But Pipeline Continued

**Cause:** Non-blocking R analysis (default behavior)

**Check logs:**
```bash
cat /path/to/species/processed_files/r_analysis_error.log
```

**Common issues:**
- Missing R packages (run install script)
- Incorrect sample names in metadata
- No GTF file (if annotation enabled)

**To require R analysis:**
```bash
sbatch --export=SPECIES_NAME=my_species,R_ANALYSIS_REQUIRED=true multi_species_atac_pipeline_v2.sh
```

---

## Advanced Configuration

### Environment Variables

Set persistent configuration in `~/.bashrc`:

```bash
# IISAGE R Analysis Configuration
export IISAGE_BASE_DIR="/path/to/atac_projects"
export IISAGE_EMAIL="your.email@institution.edu"
export R_MODULE="R/4.3.2"
export R_SKIP_MOTIFS=false
```

### Custom R Script Location

```bash
export R_SCRIPT_PATH="/custom/location/generate_atac_figs_v2.R"
sbatch --export=SPECIES_NAME=my_species multi_species_atac_pipeline_v2.sh
```

---

## Performance Tips

### 1. Skip Motifs for Speed

Motif analysis (HOMER) is the slowest step:

```bash
# Save ~30-60 minutes per sample
sbatch --export=SPECIES_NAME=my_species,R_SKIP_MOTIFS=true multi_species_atac_pipeline_v2.sh
```

### 2. Use Resume Mode

If pipeline interrupted, checkpoint system skips completed stages:

```bash
# Automatically resumes from last successful stage
sbatch --export=SPECIES_NAME=my_species multi_species_atac_pipeline_v2.sh
```

### 3. Array Jobs for Large Projects

Process samples in parallel:

```bash
# 10 samples in parallel (10x faster than sequential)
sbatch --array=1-10 --export=SPECIES_NAME=my_species multi_species_atac_pipeline_v2.sh
```

R analysis runs automatically in main job after all array tasks complete.

---

## Integration with Downstream Analysis

### Using Output Files

```r
# Load peak annotation results
library(readr)
annot <- read_csv("figures/sample1_peak_annotation.csv")

# Load peak features summary
features <- read_csv("figures/peak_features_summary.csv")

# Load consensus peaks in R
library(rtracklayer)
consensus <- import("figures/consensus_peaks_all.bed")
```

### Further Analysis

After R integration completes, you can:

1. **Customize visualizations:** Modify `generate_atac_figs_v2.R` and rerun
2. **Differential analysis:** Use DiffBind with generated data
3. **Pathway enrichment:** Use genes from peak annotations
4. **Custom motif analysis:** Use consensus peaks with different tools

---

## FAQ

**Q: Do I need to install R packages manually?**

A: No, `generate_atac_figs_v2.R` has automatic installation at the top.

**Q: Can I run R analysis on multiple species at once?**

A: Yes! R script processes all species in BASE_DIR:

```bash
Rscript generate_atac_figs_v2.R --base_dir=/path/to/base
```

**Q: What if some samples fail peak calling?**

A: R analysis runs on successfully processed samples. Check metadata matches available peak files.

**Q: How much disk space does R analysis use?**

A: Minimal (~50-100MB per species). Figures are PNG format, tables are CSV.

**Q: Can I disable email notifications for R analysis only?**

A: R analysis uses main `SEND_EMAIL_UPDATES` setting. To disable all emails:

```bash
sbatch --export=SPECIES_NAME=my_species,SEND_EMAIL_UPDATES=false multi_species_atac_pipeline_v2.sh
```

---

## Example Workflow

Complete analysis from raw data to figures:

```bash
# 1. Set environment
export IISAGE_BASE_DIR="/mnt/scratch/$USER/atac_projects"
export IISAGE_EMAIL="$USER@msu.edu"

# 2. Create sample metadata
cat > $IISAGE_BASE_DIR/my_species/sample_metadata.csv << EOF
sampID,condition,replicate
ctrl_1,control,1
ctrl_2,control,2
treat_1,treatment,1
treat_2,treatment,2
EOF

# 3. Run complete pipeline (with array jobs for speed)
cd /path/to/iisage_atac-seq/scripts
sbatch --array=1-4 --export=SPECIES_NAME=my_species multi_species_atac_pipeline_v2.sh

# 4. Monitor progress
squeue -u $USER

# 5. Check results when complete
ls -lh $IISAGE_BASE_DIR/my_species/figures/
firefox $IISAGE_BASE_DIR/my_species/figures/*.png
```

**Expected runtime:**
- With array jobs (4 samples): ~3-4 hours
- Sequential: ~10-12 hours
- R analysis: ~15-30 minutes (depends on motif analysis)

---

For more details, see:
- Main pipeline documentation: `docs/v2_upgrade_guide.md`
- R analysis details: `docs/r_analysis_guide.md`
- Troubleshooting: `docs/troubleshooting.md`

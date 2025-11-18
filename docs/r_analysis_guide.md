# IISAGE ATAC-seq R Analysis Guide

## Overview

The R analysis scripts provide downstream analysis of ATAC-seq peak data, including visualization, annotation, and motif discovery.

---

## Version Comparison

| Feature | v1.0 | v2.0 |
|---------|------|------|
| Venn diagrams | ✓ | ✓ |
| Correlation plots | ✓ | ✓ |
| Heatmaps | ✓ | ✓ |
| Peak annotation | ✗ | ✓ |
| Genomic distribution | ✗ | ✓ |
| Distance to TSS | ✗ | ✓ |
| Motif analysis (HOMER) | ✗ | ✓ |
| Peak feature tables | ✗ | ✓ |
| Consensus peak sets | ✗ | ✓ |
| Summary statistics | ✗ | ✓ |
| Configurable parameters | ✗ | ✓ |
| Command-line interface | ✗ | ✓ |

---

## Installation

### Required R Packages

```r
# Install Bioconductor if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Core packages (required)
BiocManager::install(c(
  "DiffBind",
  "GenomicRanges",
  "GenomicFeatures"
))

install.packages(c(
  "pheatmap",
  "readr",
  "dplyr",
  "ggplot2"
))

# Annotation packages (highly recommended)
BiocManager::install(c(
  "ChIPseeker",
  "rtracklayer"
))

# Motif analysis packages (optional)
BiocManager::install("memes")
```

### HOMER Installation (Optional but Recommended)

For motif analysis, install HOMER:

```bash
# Download and install HOMER
cd ~/software
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install

# Add to PATH in ~/.bashrc
export PATH=$PATH:~/software/homer/bin

# Install organism-specific data (example for human)
perl ~/software/homer/configureHomer.pl -install hg38
```

---

## Configuration Requirements

### 1. Sample Metadata File

Each species directory must contain `sample_metadata.csv`:

```csv
sampID,condition,replicate
sample1,control,1
sample2,control,2
sample3,treatment,1
sample4,treatment,2
```

**Column descriptions:**
- `sampID`: Sample identifier (must match peak file names)
- `condition`: Experimental condition
- `replicate`: Biological replicate number

### 2. Config.yaml Additions (Optional)

Add these to your `config.yaml` for enhanced analysis:

```yaml
# Existing required fields
SPECIES_NAME: my_species
GENOME_SIZE: 2.7e9
REF_GENOME_FASTA: reference/genome.fna

# NEW: For peak annotation
GTF_FILE: reference/annotation.gtf
# Or
GFF_FILE: reference/annotation.gff3

# NEW: For motif analysis (usually same as REF_GENOME_FASTA)
REF_GENOME_FASTA: reference/genome.fna
```

### 3. Directory Structure

Expected structure after running pipeline:

```
species_name/
├── config.yaml
├── sample_metadata.csv
├── reference/
│   ├── genome.fna
│   ├── annotation.gtf  (optional, for annotation)
│   └── genome_index/
├── samples/
│   ├── sample1/
│   └── sample2/
├── peaks_output/
│   ├── sample1_peaks.narrowPeak
│   ├── sample1_peaks.xls
│   ├── sample1_size_selected_100.bam
│   ├── sample2_peaks.narrowPeak
│   ├── sample2_peaks.xls
│   └── sample2_size_selected_100.bam
└── figures/  (created by R script)
```

---

## Usage

### Basic Usage (v2.0)

```bash
# From command line
cd /path/to/iisage_atac-seq/scripts
Rscript generate_atac_figs_v2.R --base_dir=/path/to/atac_seq_projects

# Skip motif analysis
Rscript generate_atac_figs_v2.R --base_dir=/path/to/atac_seq_projects --no-motifs

# Skip annotation
Rscript generate_atac_figs_v2.R --base_dir=/path/to/atac_seq_projects --no-annotation
```

### From R Console

```r
# Load the script
source("generate_atac_figs_v2.R")

# Run complete analysis
generate_comprehensive_analysis("/path/to/atac_seq_projects")

# Run without motif analysis
generate_comprehensive_analysis("/path/to/atac_seq_projects",
                                run_motifs = FALSE)

# Run without annotation
generate_comprehensive_analysis("/path/to/atac_seq_projects",
                                run_annotation = FALSE)
```

### Legacy Usage (v1.0)

```r
source("generate_atac_figs.R")
generate_figures_modular("atac_sec/")
```

---

## Output Files

The script creates a `figures/` directory in each species folder with the following outputs:

### 1. Core Visualizations (v1.0 + v2.0)

| File | Description |
|------|-------------|
| `venn_all.png` | Venn diagram showing peak overlaps across all samples |
| `venn_<condition>.png` | Condition-specific Venn diagrams |
| `correlation.png` | Sample correlation plot based on peak overlap |
| `heatmap_base.png` | Basic heatmap of peak intensities |
| `heatmap_pheatmap.png` | Enhanced clustered heatmap |

### 2. Peak Annotation Files (v2.0 NEW)

| File | Description |
|------|-------------|
| `<sample>_peak_annotation.csv` | Complete peak annotation table with genomic features |
| `<sample>_annotation_pie.png` | Pie chart of genomic feature distribution |
| `<sample>_annotation_bar.png` | Bar chart of genomic feature distribution |
| `<sample>_distance_to_tss.png` | Distribution of peak distances to TSS |

### 3. Summary Tables (v2.0 NEW)

| File | Description |
|------|-------------|
| `peak_features_summary.csv` | Comprehensive peak statistics per sample |
| `genomic_distribution_summary.csv` | Percentage of peaks in each genomic region |
| `peak_counts_by_sample.png` | Bar plot of total peaks per sample |
| `peak_width_distribution.png` | Peak width statistics visualization |

### 4. Comparative Analysis (v2.0 NEW)

| File | Description |
|------|-------------|
| `genomic_distribution_comparison.png` | Side-by-side genomic distribution across samples |
| `distance_to_tss_comparison.png` | Comparative TSS distance distribution |

### 5. Consensus Peak Sets (v2.0 NEW)

| File | Description |
|------|-------------|
| `consensus_peaks_all.bed` | Union of all peaks across samples (BED format) |
| `consensus_peaks_<condition>.bed` | Consensus peaks per condition |

### 6. Motif Analysis (v2.0 NEW)

| Directory | Description |
|-----------|-------------|
| `<sample>_motifs/` | HOMER motif analysis results |
| `<sample>_motifs/homerResults.html` | Interactive HTML report with discovered motifs |
| `<sample>_motifs/knownResults.html` | Known motif enrichment results |

---

## Understanding the Outputs

### Peak Annotation CSV Format

Example columns in `*_peak_annotation.csv`:

```csv
seqnames,start,end,width,strand,annotation,distanceToTSS,geneId,transcriptId
chr1,1000000,1000500,500,*,Promoter (2-3kb),-2500,GENE001,TRANS001
chr1,2000000,2000300,300,*,Intron,15000,GENE002,TRANS002
```

**Key columns:**
- `annotation`: Genomic feature type (Promoter, Exon, Intron, Intergenic, etc.)
- `distanceToTSS`: Distance to nearest transcription start site (bp)
- `geneId`: Associated gene identifier
- `transcriptId`: Associated transcript identifier

### Peak Features Summary CSV

Example columns in `peak_features_summary.csv`:

```csv
sample,condition,replicate,total_peaks,mean_peak_width,median_peak_width,mean_signal,median_signal,mean_qvalue,peaks_high_confidence,peaks_strong_signal
sample1,control,1,25000,350,300,45.2,38.5,25.3,15000,6250
sample2,control,2,28000,340,295,48.1,40.2,27.1,17000,7000
```

**Interpretation:**
- `total_peaks`: Number of MACS2-called peaks
- `mean/median_peak_width`: Average peak size in bp
- `mean/median_signal`: Average signal intensity
- `mean_qvalue`: Average -log10(q-value)
- `peaks_high_confidence`: Peaks with q-value > 10
- `peaks_strong_signal`: Peaks in top 25% by signal

### Genomic Distribution Summary

Example columns in `genomic_distribution_summary.csv`:

```csv
sample,promoter,exon,intron,intergenic,mean_dist_to_tss
sample1,35.2,8.5,28.3,28.0,12500
sample2,38.1,9.2,26.7,26.0,11800
```

**Interpretation:**
- Values are percentages (sum to 100%)
- `promoter`: % of peaks in promoter regions (±3kb TSS by default)
- `exon/intron/intergenic`: % of peaks in respective features
- `mean_dist_to_tss`: Average distance to nearest TSS (bp)

### HOMER Motif Results

Check `<sample>_motifs/homerResults.html` for:
- **De novo motifs**: Novel motifs discovered in your peaks
- **Motif sequences**: Consensus sequences (e.g., ATGCAAAT)
- **P-values**: Statistical significance of enrichment
- **% targets**: Percentage of peaks containing motif
- **Comparison to known motifs**: Similar transcription factor binding sites

---

## Customization

### Adjusting TSS Regions

Edit the script constants:

```r
# Default: ±3kb around TSS
DEFAULT_PROMOTER_UPSTREAM <- 3000
DEFAULT_PROMOTER_DOWNSTREAM <- 3000
DEFAULT_TSS_REGION <- c(-3000, 3000)

# Change to ±5kb:
DEFAULT_PROMOTER_UPSTREAM <- 5000
DEFAULT_PROMOTER_DOWNSTREAM <- 5000
DEFAULT_TSS_REGION <- c(-5000, 5000)
```

### Adjusting Figure Quality

```r
# Default DPI and sizes
FIG_DPI <- 300          # Change to 600 for publication quality
FIG_WIDTH <- 7          # Inches
FIG_HEIGHT <- 7
FIG_WIDTH_WIDE <- 10    # For multi-panel plots
```

### HOMER Motif Parameters

Edit the HOMER command in `run_homer_motif_analysis()`:

```r
# Current: 200bp regions centered on peak summit
homer_cmd <- sprintf(
  "findMotifsGenome.pl '%s' '%s' '%s' -size 200 -mask",
  peaks_file,
  genome_fasta,
  motif_dir
)

# Use peak size instead of fixed 200bp:
homer_cmd <- sprintf(
  "findMotifsGenome.pl '%s' '%s' '%s' -size given -mask",
  peaks_file,
  genome_fasta,
  motif_dir
)

# Add motif length constraints:
homer_cmd <- sprintf(
  "findMotifsGenome.pl '%s' '%s' '%s' -size 200 -mask -len 8,10,12",
  peaks_file,
  genome_fasta,
  motif_dir
)
```

---

## Troubleshooting

### Error: "ChIPseeker not available"

**Solution:**
```r
BiocManager::install("ChIPseeker")
```

If that fails, try:
```r
BiocManager::install("ChIPseeker", force = TRUE)
```

### Error: "No GTF/GFF annotation file found"

**Solution:**

1. Download annotation for your organism:
   - Ensembl: https://www.ensembl.org/
   - NCBI: https://www.ncbi.nlm.nih.gov/genome/
   - UCSC: https://genome.ucsc.edu/

2. Add to `config.yaml`:
   ```yaml
   GTF_FILE: reference/annotation.gtf
   ```

3. Or place in `reference/` directory (auto-detected)

### Error: "HOMER not found"

**Solution:**

Either:
1. Install HOMER (see Installation section)
2. Skip motif analysis: `--no-motifs` flag

### Warning: "no valid samples, skipping"

**Causes:**
- Missing `sample_metadata.csv`
- Peak files don't match `sampID` in metadata
- BAM files missing

**Solution:**

Check file naming:
```bash
# Metadata says "sample1"
# These files must exist:
peaks_output/sample1_peaks.narrowPeak
peaks_output/sample1_peaks.xls
peaks_output/sample1_size_selected_100.bam
```

### Error: Package loading issues

**Solution:**

Test package availability:
```r
# Check if packages load
library(DiffBind)
library(ChIPseeker)
library(GenomicFeatures)

# If error, reinstall
BiocManager::install(c("DiffBind", "ChIPseeker", "GenomicFeatures"))
```

---

## Best Practices

### 1. Quality Control

Before running analysis, check:
- All samples have similar peak counts (±50%)
- Peak files are not empty
- BAM files are indexed (`.bai` files present)

### 2. Biological Replicates

- Use at least 2 replicates per condition
- More replicates = better statistical power
- Label correctly in `sample_metadata.csv`

### 3. Genome Annotation

- Use annotation matching your reference genome version
- GTF preferred over GFF for compatibility
- Compressed files (`.gz`) are supported

### 4. Motif Analysis

- HOMER requires 50+ peaks for reliable results
- Larger peak sets (1000+) give better motif discovery
- Use appropriate genome background (same as alignment)

### 5. Reproducibility

Save session info:
```r
# After running analysis
sink("R_session_info.txt")
sessionInfo()
sink()
```

---

## Integration with Pipeline

### Automated Integration (Future)

The v2.0 bash pipeline will eventually call this R script automatically. For now, run manually after pipeline completion:

```bash
# 1. Run pipeline
sbatch --export=SPECIES_NAME=my_species multi_species_atac_pipeline_v2.sh

# 2. Wait for completion, then run R analysis
Rscript generate_atac_figs_v2.R --base_dir=$IISAGE_BASE_DIR
```

### Batch Processing

Process multiple species:

```bash
# All species at once (if sample_metadata.csv exists for each)
Rscript generate_atac_figs_v2.R --base_dir=/path/to/atac_seq_projects

# Single species
Rscript generate_atac_figs_v2.R --base_dir=/path/to/atac_seq_projects/species1
```

---

## Example Workflow

Complete analysis workflow:

```bash
# 1. Setup environment
export IISAGE_BASE_DIR="/mnt/scratch/$USER/atac_seq_projects"
export IISAGE_EMAIL="$USER@msu.edu"

# 2. Run ATAC-seq pipeline
cd /path/to/iisage_atac-seq/scripts
sbatch --array=1-10 --export=SPECIES_NAME=my_species multi_species_atac_pipeline_v2.sh

# 3. Monitor job
squeue -u $USER

# 4. After completion, create sample_metadata.csv
cat > $IISAGE_BASE_DIR/my_species/sample_metadata.csv << EOF
sampID,condition,replicate
sample1,control,1
sample2,control,2
sample3,treatment,1
sample4,treatment,2
EOF

# 5. Run R analysis
Rscript generate_atac_figs_v2.R --base_dir=$IISAGE_BASE_DIR

# 6. Review results
cd $IISAGE_BASE_DIR/my_species/figures
ls -lh
firefox sample1_motifs/homerResults.html  # View motif results
```

---

## Output Interpretation Guide

### What to Look For

**1. Venn Diagrams:**
- High overlap (>50%) = good reproducibility between replicates
- Low overlap (<30%) = check sample quality or biological variability

**2. Correlation Heatmap:**
- Replicates should cluster together
- Pearson r > 0.8 = good correlation
- Check for batch effects (samples clustering by prep date, not biology)

**3. Peak Features:**
- Typical ATAC-seq: 20,000-100,000 peaks
- Mean width: 200-500 bp
- Compare similar tissues/conditions for consistency

**4. Genomic Distribution:**
- ATAC-seq typically enriched in promoters (30-50%)
- High intergenic % might indicate enhancers
- Compare to published datasets for your organism/tissue

**5. Distance to TSS:**
- ATAC-seq peaks cluster around TSS (±2kb)
- Bimodal distribution normal (promoters + distal regulatory)
- Compare conditions to identify differential accessibility

**6. Motifs:**
- Top motifs should match expected transcription factors
- Check literature for your tissue/condition
- Novel motifs warrant further investigation

---

## Advanced Usage

### Custom Annotation Databases

Use organism-specific databases:

```r
# For mouse
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# For human
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Modify script to use specific TxDb instead of loading from GTF
```

### Gene Ontology Enrichment

Extend peak annotation with GO analysis:

```r
# After annotating peaks
library(clusterProfiler)

# Extract genes near peaks
peak_genes <- unique(anno_df$geneId)

# GO enrichment
ego <- enrichGO(gene = peak_genes,
                OrgDb = org.Hs.eg.db,  # Change for organism
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# Plot
dotplot(ego, showCategory=20)
```

---

## Citation

If you use this pipeline in your research, please cite:

- **DiffBind**: Ross-Innes CS, et al. (2012) Genome Research
- **ChIPseeker**: Yu G, et al. (2015) Bioinformatics
- **HOMER**: Heinz S, et al. (2010) Molecular Cell
- **MACS2**: Zhang Y, et al. (2008) Genome Biology

And acknowledge:
> "ATAC-seq analysis was performed using the IISAGE pipeline (https://github.com/...)"

---

*For questions or issues, contact: Louis Paul Decena-Segarra, PhD*

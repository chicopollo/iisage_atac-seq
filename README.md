# IISAGE: Multi-Species ATAC-seq Pipeline

Comprehensive and standardized ATAC-seq analysis and figure generation pipeline for SLURM-based HPCC environments.
(so far this has only been tested in MSU's HPCC)

**Current Version:** 2.5.0

## Overview

The IISAGE ATAC-seq pipeline provides end-to-end analysis of chromatin accessibility data with:
- **Two-job workflow**: Parallel sample processing + integrated downstream analysis
- **Smart resume capabilities**: Automatic detection and reuse of intermediate files
- **Multi-species support**: Process multiple species with species-specific configurations
- **Comprehensive outputs**: From raw reads to publication-ready figures

## Features

### Pipeline Capabilities
- ✓ **Parallel processing**: SLURM array jobs for 16x faster sample processing
- ✓ **Automated workflow**: QC → trimming → alignment → peak calling → visualization → R analysis
- ✓ **Smart resume**: Checkpoints at every stage, validates file integrity
- ✓ **Modular stages**: Enable/disable components independently
- ✓ **R integration**: DiffBind, deepTools, peak annotation, motif analysis

### Performance Optimizations (v2.5)
- **Parallelized peak calling**: 16 samples processed simultaneously (16x speedup)
- **Optimized two-job workflow**:
  - Job 1 (Array): Per-sample processing in parallel
  - Job 2 (Analysis): Multi-sample comparisons after array completion
- **Comprehensive resume logic**: Skips completed steps (BAM, SAM, trimmed FASTQ, merged FASTQ, etc)
- **File validation**: Detects and reprocesses corrupted/incomplete files

### Module Conflict Resolution
- **Stage-specific module loading**: Separate modules are called at required times

## Quick Start

### 1. Clone Repository
```bash
git clone https://github.com/iisage/iisage_atac-seq.git
cd iisage_atac-seq
```

### 2. Set Up R Environment (One-time)

Install required R packages:
```bash
cd scripts

# Option A: Interactive installation (recommended)
module load R/4.2.2-foss-2022b
R
source("install_r_packages.R")

# Option B: Submit as SLURM job
sbatch install_r_packages_job.sh
```

**Note**: R package installiation should be done into user specific environment

### 3. Organize Your Data

```
/path/to/atac_seq_projects/
├── Species_name/
│   ├── config.yaml                    # Species-specific configuration
│   ├── sample_metadata.csv            # Experimental design (CSV format)
│   ├── samples/
│   │   ├── sample1/
│   │   │   ├── *_1.fq.gz             # Forward reads
│   │   │   └── *_2.fq.gz             # Reverse reads
│   │   ├── sample2/
│   │   └── sample3/
│   └── reference/
│       └── genome.fasta               # Reference genome
```

### 4. Configure Your Species

**config.yaml** (place in species directory):
```yaml
SPECIES_NAME: species_name
GENOME_SIZE: 2.4e9
REF_GENOME_FASTA: reference/genome.fasta
ADAPTERS_PATH: /path/to/adapters.fa
TRIM_SETTINGS: "ILLUMINACLIP:/path/to/adapters.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:25"
```
**Note**: adapters.fa is a multifasta file with adapter sequence, provided by sequencing core.

**sample_metadata.csv** (required for R analysis):
```csv
sampID,TrueID,condition,sex,age_category,replicate
sample1,ID1,treatment,M,Y,1
sample2,ID2,treatment,F,Y,2
sample3,ID3,control,M,O,1
```

**Required columns**: `sampID`, `condition`, `replicate`
**Optional columns**: `sex`, `age_category`, `TrueID`, custom fields

`sampID`: Sequencing core sample identification ID
`condition`: Indicates what kind of sample it is (i.e: I use it to indicate cell tissue type)

More custom fields can be aded to the metadata file.

### 5. Deploy Pipeline

Update paths in `atac_seq_pipeline_deployment.sh`:
```bash
# Edit line 27-28
BASE_DIR=/mnt/scratch/your_user/atac_seq_projects"
PIPELINE_SCRIPT="$SCRIPT_DIR/multi_species_atac_pipeline_v2.5.sh"
```

**Validate configuration**:
```bash
./atac_seq_pipeline_deployment.sh --validate --species Species_name
```

**Deploy single species**:
```bash
./atac_seq_pipeline_deployment.sh --species Species_name
```

**Deploy all valid species**:
```bash
./atac_seq_pipeline_deployment.sh --deploy-all
```

### 6. Monitor Jobs

```bash
# Check job status
squeue -u $USER

# View array job logs (16 tasks)
tail -f atac-seq_pipeline-JOBID_*.SLURMout

# Check detailed job info
sacct -j JOBID,JOBID2 --format=JobID,JobName,State,ExitCode,Elapsed
```

## Pipeline Stages

### Job 1: Array Processing (Parallel)
1. **Trimming**: Adapter removal with Trimmomatic
2. **Alignment**: Bowtie2 alignment to reference genome
3. **BAM Processing**: SAM to BAM conversion, sorting, indexing
4. **Size Selection**: Fragment size filtering (<100bp removed)
5. **Peak Calling**: MACS2 peak calling (per sample, in parallel)

### Job 2: Analysis (Sequential)
6. **Peak Verification**: Confirm all samples have peaks
7. **DeepTools Visualization**:
   - Correlation heatmaps
   - PCA plots
   - Coverage tracks (bigWig)
8. **R Analysis** (DiffBind):
   - Differential accessibility
   - Venn diagrams
   - Peak annotation
   - Motif enrichment (if HOMER available)

## Requirements

### SLURM Modules
Verify modules are available and note exact versions:
```bash
module spider Trimmomatic
module spider Bowtie2
module spider SAMtools
module spider MACS2
module spider deepTools
module spider R
```

**Recommended versions** (MSU HPCC):
- Trimmomatic (any recent version)
- Bowtie2/2.5.1-GCC-12.2.0
- SAMtools/1.17-GCC-12.2.0
- MACS2/2.2.9.1-foss-2022b
- deepTools/3.5.5-foss-2023a
- R/4.2.2-foss-2022b

### R Packages (Auto-installed)

**Required**:
- DiffBind, GenomicRanges, GenomicFeatures, rtracklayer, Rsamtools
- pheatmap, ggplot2, dplyr, readr

**Optional**:
- ChIPseeker (peak annotation)
- ATACseqQC (ATAC-seq specific QC)
- clusterProfiler (pathway enrichment)

### System Requirements

**Per array task**:
- CPUs: 4 (default)
- Memory: 16GB (4GB/CPU)
- Time: ~1-2 hours per sample

**Analysis job**:
- CPUs: 8
- Memory: 64GB (8GB/CPU)
- Time: ~2-4 hours (depends on sample count)

## Output Directory Structure

```
Species_name/
├── processed_files/
│   ├── *_sorted.bam                   # Aligned BAM files
│   └── *_sorted.bam.bai               # BAM indices
├── peaks_output/
│   ├── *_size_selected_100.bam        # Size-filtered BAMs
│   ├── *_peaks.narrowPeak             # MACS2 peaks
│   ├── *_summits.bed                  # Peak summits
│   └── consensus_peaks.bed            # Merged consensus peaks
├── deeptools_output/
│   ├── correlation_heatmap.png        # Sample correlations
│   ├── pca_plot.png                   # PCA analysis
│   └── *.bigWig                       # Coverage tracks
├── figures/
│   ├── venn_diagrams/                 # Peak overlaps
│   ├── diffbind_results/              # Differential accessibility
│   ├── peak_annotation/               # Genomic distribution
│   └── motif_analysis/                # Enriched motifs (if available)
└── checkpoints/
    └── *.done                         # Resume checkpoints
```

## Resume Capabilities

The pipeline automatically resumes from the most processed state:

1. **Valid BAM exists** → Skip all preprocessing
2. **Valid SAM exists** → Resume from BAM conversion
3. **Trimmed FASTQs exist** → Resume from alignment
4. **Merged FASTQs exist** → Resume from trimming
5. **Size-selected BAM exists** → Skip size selection
6. **Peaks exist** → Skip peak calling

**File validation**: Corrupted/incomplete files are detected and reprocessed automatically.

## Test Dataset

Request test data from: abroniko@msu.edu

**Species included**:
- *Chrysemys picta* (painted turtle)
- *Thamnophis elegans* (garter snake)
- *Hemidactylus turcicus* (Mediterranean house gecko)

**For Bronikowski Lab members** (Carson access):
```bash
# From HPCC
rsync -r $USER@carson.kbs.msu.edu:/data/grpdata/broniko_lab/IISAGE_atac_seq/test_data .
```

## Troubleshooting

### Module Conflicts

**Problem**: `foss-2022b` and `foss-2023a` incompatible

**Solution**: Pipeline uses stage-specific module loading:
- Stage 1: Trimmomatic, Bowtie2, SAMtools (foss-2022b compatible)
- Stage 2: MACS2 (requires foss-2022b)
- Stage 3: deepTools (requires foss-2023a)

### Job Failures

**Check logs**:
```bash
# Array task errors
grep ERROR atac-seq_pipeline-JOBID_*.SLURMout

# Analysis job errors
grep ERROR atac-seq_pipeline-JOBID.SLURMout

# SLURM errors (memory, timeout)
cat atac-seq_pipeline-JOBID.SLURMerr
```

**Common issues**:
1. **Out of memory**: Increase `--mem-per-cpu` in deployment script
2. **Timeout**: Increase time limit for array jobs
3. **Module not found**: Check module names/versions with `module spider`
4. **Corrupted SAM/BAM**: Delete and re-run (pipeline auto-detects)

### R Package Installation Fails

**Issue**: BiocManager packages fail to install

**Solutions**:
1. Install problematic packages manually from live session:
   ```R
   BiocManager::install("rtracklayer", force=TRUE, type="source")
   ```

### Deployment Script Issues

**Validate first**:
```bash
./atac_seq_pipeline_deployment.sh --validate --species Species_name
```

**Check**:
- `config.yaml` exists and has required fields
- `sample_metadata.csv` is comma-separated (not tab-separated)
- Sample directories contain paired FASTQ files (`*_1.fq.gz`, `*_2.fq.gz`)
- Reference genome path is correct

## Version History

### v2.5.0 (Current)
- Peak calling parallelized in array stage (16x speedup)
- SAM file validation before resume
- Smart resume from intermediate FASTQ files
- Comprehensive checkpoint system

### v2.4.0
- Size selection moved to array stage
- Two-job workflow implementation
- Optimized resource allocation

### v2.3.0
- Smart array job coordination for index building
- File locking mechanism

### v2.2.0
- Split module loading to avoid conflicts
- Stage-specific module management

### v2.0.0
- YAML-based configuration
- SLURM array job support
- Integrated R analysis

## Citation

If you use this pipeline, please cite:

> Decena-Segarra, L.P., Bronikowski, A.M. (2025). IISAGE: Integrated Multi-Species ATAC-seq Analysis Pipeline. GitHub repository. https://github.com/iisage/iisage_atac-seq

## Support

- **Issues**: https://github.com/yourusername/iisage_atac-seq/issues
- **Contact**: decenalo@msu.edu

## License

MIT License - See LICENSE file for details

---

**Author**: Louis Paul Decena-Segarra, PhD
**Lab**: Bronikowski Lab, Michigan State University
**Last Updated**: 2025-11-18

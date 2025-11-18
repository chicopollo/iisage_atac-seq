# IISAGE ATAC-seq Pipeline v2.0 - Upgrade Guide

## What's New in v2.0

Version 2.0 introduces significant improvements to robustness, configurability, and performance.

### Major Features

1. **Environment-based Configuration** - No more hard-coded paths
2. **Safer YAML Parsing** - Removed `eval` security vulnerability
3. **Module Validation** - Ensures all required software is available before running
4. **SLURM Array Jobs** - Parallel sample processing for massive speedups
5. **Pipeline Resume** - Skip completed stages automatically
6. **Disk Space Checks** - Prevents failures from running out of space
7. **Comprehensive Error Checking** - Better error messages and validation
8. **Configurable Parameters** - All magic numbers now configurable
9. **Dry-Run Mode** - Preview what will happen without executing
10. **Version Tracking** - Records all software versions for reproducibility

---

## Configuration Changes

### Environment Variables

Instead of editing the script, set these environment variables:

```bash
# Required: Base directory for all species data
export IISAGE_BASE_DIR="/path/to/your/atac_seq_projects"

# Optional: Working directory for SLURM jobs (default: /mnt/scratch/$USER)
export IISAGE_WORK_DIR="/path/to/scratch"

# Optional: Email for notifications (default: $USER@msu.edu)
export IISAGE_EMAIL="your.email@institution.edu"
```

Add these to your `~/.bashrc` or `~/.bash_profile` for persistence.

### YAML Configuration Enhancements

New optional parameters you can add to `config.yaml`:

```yaml
# Existing required parameters
SPECIES_NAME: species_name
GENOME_SIZE: 2.7e9
REF_GENOME_FASTA: reference/genome.fna

# NEW: DeepTools parameters
TSS_UPSTREAM: 2000          # bp upstream of TSS (default: 2000)
TSS_DOWNSTREAM: 2000        # bp downstream of TSS (default: 2000)
BIGWIG_BINSIZE: 10          # BigWig bin size in bp (default: 10)
BIGWIG_NORMALIZE: RPKM      # Normalization method: RPKM, CPM, BPM, RPGC (default: RPKM)

# Existing optional parameters (now with documented defaults)
THREADS: 16                 # CPU threads (default: 16)
MAX_INSERT_SIZE: 2000       # Maximum fragment size (default: 2000)
SIZE_SELECTION_CUTOFF: 100  # Nucleosome-free cutoff (default: 100)
MIN_MAPQ: 5                 # Minimum mapping quality (default: 5)
```

---

## Usage Examples

### 1. Basic Usage (No Changes)

```bash
# Single species
sbatch --export=SPECIES_NAME=cpicta multi_species_atac_pipeline_v2.sh

# All species
sbatch --export=DEPLOY_ALL=true multi_species_atac_pipeline_v2.sh
```

### 2. Dry-Run Mode (NEW)

Preview what the pipeline will do without executing:

```bash
sbatch --export=SPECIES_NAME=cpicta,DRY_RUN=true multi_species_atac_pipeline_v2.sh
```

Check the output file to see all commands that would be executed.

### 3. Array Job Mode (NEW - RECOMMENDED)

For **significantly faster** processing when you have multiple samples:

**Step 1:** Count your samples
```bash
ls -d /path/to/species/samples/*/ | wc -l
# Example output: 12
```

**Step 2:** Submit array job
```bash
# If you have 12 samples, use --array=1-12
sbatch --array=1-12 --export=SPECIES_NAME=cpicta multi_species_atac_pipeline_v2.sh
```

This will:
- Process all 12 samples in **parallel** (up to SLURM's array limit)
- Each sample gets its own compute node
- 12x faster than sequential processing
- After all samples complete, run peak calling and visualization

### 4. Resume After Failure (NEW)

If the pipeline fails or is cancelled, just resubmit:

```bash
sbatch --export=SPECIES_NAME=cpicta multi_species_atac_pipeline_v2.sh
```

The pipeline will:
- Skip already-completed stages
- Resume from the last checkpoint
- Save time and compute resources

To **force** re-running everything:
```bash
# Remove checkpoint directory
rm -rf /path/to/species/.checkpoints
sbatch --export=SPECIES_NAME=cpicta multi_species_atac_pipeline_v2.sh
```

### 5. Custom Configuration

Override defaults without editing YAML:

```bash
sbatch --export=SPECIES_NAME=cpicta,IISAGE_BASE_DIR=/custom/path,IISAGE_EMAIL=me@email.com \
       multi_species_atac_pipeline_v2.sh
```

---

## New Output Files

### Software Versions File

`processed_files/software_versions.txt` now contains:
- Pipeline version
- All tool versions (Trimmomatic, Bowtie2, SAMtools, MACS2, deepTools)
- SLURM job information
- Configuration parameters
- Timestamp and host information

Example:
```
=== Software Versions ===
Pipeline: IISAGE ATAC-seq Pipeline v2.0.0
Date: 2025-11-14 10:30:45
Host: dev-intel18-k80
User: decenalo

--- Core Tools ---
Trimmomatic: 0.39
Bowtie2: version 2.4.5
SAMtools: 1.16
MACS2: 2.2.9.1
deepTools: bamCoverage 3.5.1

--- Environment ---
SLURM_JOB_ID: 12345678
CPUs: 16
Genome: /path/to/genome.fna
Genome Size: 2.7e9
```

### Checkpoint Files

`.checkpoints/` directory tracks completed stages:
- `bowtie2_index.done`
- `sample_<name>_aligned.done`
- `sample_<name>_peaks.done`
- `deeptools_complete.done`

### Enhanced Summary File

`processed_files/pipeline_summary.txt` now includes:
- Pipeline version
- All configuration parameters
- Software versions reference
- Detailed statistics

---

## Migration Guide

### From v1.x to v2.0

1. **Set environment variables** (one-time setup):
   ```bash
   echo 'export IISAGE_BASE_DIR="/mnt/scratch/decenalo/atac_seq_projects"' >> ~/.bashrc
   echo 'export IISAGE_EMAIL="decenalo@msu.edu"' >> ~/.bashrc
   source ~/.bashrc
   ```

2. **Update SLURM submission scripts**:

   Old:
   ```bash
   sbatch --export=SPECIES_NAME=cpicta multi_species_atac_pipeline.sh
   ```

   New (same, but optionally use arrays):
   ```bash
   # Sequential (same as before)
   sbatch --export=SPECIES_NAME=cpicta multi_species_atac_pipeline_v2.sh

   # Array mode (recommended for speed)
   sbatch --array=1-N --export=SPECIES_NAME=cpicta multi_species_atac_pipeline_v2.sh
   ```

3. **Optional: Update config.yaml** with new parameters:
   ```yaml
   TSS_UPSTREAM: 3000        # If you want larger TSS regions
   BIGWIG_BINSIZE: 20        # If you want coarser resolution
   ```

4. **Test with dry-run**:
   ```bash
   sbatch --export=SPECIES_NAME=cpicta,DRY_RUN=true multi_species_atac_pipeline_v2.sh
   ```

---

## Performance Improvements

### Array Jobs vs Sequential

Example with 10 samples, each taking ~2 hours:

| Mode | Time | Notes |
|------|------|-------|
| **v1.0 Sequential** | 20 hours | All samples processed one-by-one |
| **v2.0 Sequential** | 20 hours | Same as v1.0, but with resume capability |
| **v2.0 Array** | ~2 hours | All samples processed in parallel |

**Speedup**: Up to **10x faster** with array jobs!

### Resume Capability

If a job fails at hour 18 out of 20:

| Mode | Re-run Time |
|------|-------------|
| **v1.0** | 20 hours (start from scratch) |
| **v2.0** | 2 hours (resume from checkpoint) |

---

## Troubleshooting

### Q: "ERROR: Required variable GENOME_SIZE not found in config.yaml"

**A:** Ensure your `config.yaml` has proper formatting:
```yaml
GENOME_SIZE: 2.7e9  # Correct
```

Not:
```yaml
genome_size: 2.7e9  # Wrong - must be UPPERCASE
```

### Q: "ERROR: Failed to load module: Trimmomatic"

**A:** Check module availability:
```bash
module spider Trimmomatic
```

If not found, edit line 344 in the script to match your HPCC's module names.

### Q: "ERROR: Insufficient disk space. Available: 50GB, Required: 100GB"

**A:** Either:
1. Free up space in your scratch directory
2. Change the requirement (edit `MIN_FREE_SPACE_GB` in script if 100GB is too conservative)
3. Move to a different directory with more space

### Q: Array job only processes one sample?

**A:** Make sure your array size matches your sample count:
```bash
# Count samples
ls -d /path/to/species/samples/*/ | wc -l

# Use that number in --array
sbatch --array=1-12 --export=SPECIES_NAME=cpicta multi_species_atac_pipeline_v2.sh
```

### Q: How do I force re-run everything?

**A:** Delete the checkpoint directory:
```bash
rm -rf /path/to/species/.checkpoints
```

---

## Security Improvements

### v1.0 YAML Parsing (Unsafe)
```bash
eval "export $line"  # Could execute malicious code from YAML
```

### v2.0 YAML Parsing (Safe)
```bash
set_var_safe() {
    # Validates variable names
    # Sanitizes values
    # No code execution
}
```

**Impact**: Prevents YAML injection attacks. Your config files are now safely sandboxed.

---

## Best Practices

### 1. Always Use Dry-Run First

Before processing expensive datasets:
```bash
sbatch --export=SPECIES_NAME=new_species,DRY_RUN=true multi_species_atac_pipeline_v2.sh
```

Review the `.SLURMout` file to ensure everything looks correct.

### 2. Use Array Jobs for Production

Sequential mode is good for testing, but use arrays for real data:
```bash
# Testing: sequential
sbatch --export=SPECIES_NAME=test_species multi_species_atac_pipeline_v2.sh

# Production: array
sbatch --array=1-N --export=SPECIES_NAME=real_species multi_species_atac_pipeline_v2.sh
```

### 3. Keep Software Versions Files

The `software_versions.txt` file is crucial for reproducibility. Include it when:
- Publishing results
- Sharing data
- Troubleshooting issues

### 4. Monitor Disk Space

For large projects, check disk usage regularly:
```bash
du -sh /path/to/species/*
```

Enable `CLEANUP_INTERMEDIATE=true` in the script to remove temporary files.

### 5. Version Control Your Configs

Keep your `config.yaml` files in git:
```bash
cd /path/to/atac_seq_projects
git add species/*/config.yaml
git commit -m "Add ATAC-seq configs for project X"
```

---

## Feature Comparison

| Feature | v1.0 | v2.0 |
|---------|------|------|
| Hard-coded paths | ✓ | ✗ |
| Environment config | ✗ | ✓ |
| Unsafe YAML parsing | ✓ | ✗ |
| Module validation | ✗ | ✓ |
| Array job support | ✗ | ✓ |
| Resume capability | ✗ | ✓ |
| Disk space checks | ✗ | ✓ |
| Error checking | Basic | Comprehensive |
| Configurable parameters | Partial | All |
| Dry-run mode | ✗ | ✓ |
| Version tracking | ✗ | ✓ |
| Input validation | ✗ | ✓ |

---

## Getting Help

1. **Check the dry-run output** - Most issues are visible here
2. **Review SLURM error files** - `atac-species-jobid.SLURMerr`
3. **Check software versions** - `processed_files/software_versions.txt`
4. **Verify checkpoints** - `ls -la .checkpoints/`
5. **Contact** - Include pipeline version from log output

---

## Changelog

### v2.0.0 (2025-11-14)

**Added:**
- Environment variable configuration (IISAGE_BASE_DIR, IISAGE_WORK_DIR, IISAGE_EMAIL)
- Safe YAML parsing without eval
- Module loading validation
- SLURM array job support for parallel processing
- Pipeline resume/checkpoint system
- Disk space pre-flight checks
- Comprehensive error checking with detailed messages
- Configurable DeepTools parameters (TSS regions, BigWig settings)
- Dry-run mode
- Software version tracking
- Input validation (genome size, file formats)
- Bowtie2 index validation (checks all 6 index files)

**Changed:**
- Improved error messages
- Enhanced logging
- Better email notification handling
- More robust file handling

**Fixed:**
- Security vulnerability in YAML parsing
- Silent module loading failures
- Race conditions in file checking
- Incomplete Bowtie2 index validation

**Removed:**
- Hard-coded paths in script body
- Hard-coded email addresses
- Unsafe eval usage

---

## Future Enhancements

Planned for v2.1:
- Integration of `generate_atac_figs.R` into main pipeline
- MultiQC report generation
- Automated quality thresholds
- Peak annotation
- Differential accessibility analysis

---

*For questions or issues, contact: Louis Paul Decena-Segarra, PhD (decenalo@msu.edu)*

#!/bin/bash --login

##### SLURM Resource Requests #####

#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G

##### SLURM Administrative Settings #####

#SBATCH --job-name=atac-seq_pipeline
#SBATCH --output=%x-%j.SLURMout
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=CONFIGURE_EMAIL_HERE
#SBATCH --error=%x-%j.SLURMerr

####################################################################################
# Script: multi_species_atac_pipeline_v2.3.sh
# Purpose: Multi-species ATAC-seq pipeline with YAML configuration (Enhanced Version)
# Author: Louis Paul Decena-Segarra, PhD
# Version: 2.3.0
# Created: 2025-06-02
# Modified: 2025-11-17
#
# Improvements in v2.0:
#   - Configurable paths via environment variables
#   - Safer YAML parsing without eval
#   - Module loading validation
#   - SLURM array-based parallel sample processing
#   - Pipeline resume/checkpoint capability
#   - Disk space pre-flight checks
#   - Comprehensive error checking
#   - Configurable magic numbers
#   - Dry-run mode
#   - Version tracking
#   - Integrated R-based downstream analysis
#
# Improvements in v2.2:
#   - Split module loading to avoid foss toolchain conflicts
#   - Stage-specific module loading (trimming/alignment, peak calling, deepTools)
#   - Improved compatibility with HPCC module systems
#
# Improvements in v2.3:
#   - Smart array job coordination for index building
#   - File lock mechanism prevents duplicate index builds
#   - Waiter tasks automatically wait for first task to complete index
#   - Stale lock detection and cleanup
#   - Configurable timeout with helpful error messages
#
# Usage:
#   Single species: sbatch --export=SPECIES_NAME=species_name multi_species_atac_pipeline_v2.3.sh
#   All species: sbatch --export=DEPLOY_ALL=true multi_species_atac_pipeline_v2.3.sh
#   Dry run: sbatch --export=SPECIES_NAME=species_name,DRY_RUN=true multi_species_atac_pipeline_v2.3.sh
#   Array mode: sbatch --array=1-N --export=SPECIES_NAME=species_name multi_species_atac_pipeline_v2.3.sh
#   Skip R analysis: sbatch --export=SPECIES_NAME=species_name,RUN_R_ANALYSIS=false multi_species_atac_pipeline_v2.3.sh
#
# Pipeline stages: 1) Quality control & trimming, 2) Alignment & BAM processing,
#                  3) Peak calling with size selection, 4) DeepTools visualization,
#                  5) R-based downstream analysis (Venn diagrams, peak annotation, motifs)
#
# Directory structure expected:
# $IISAGE_BASE_DIR/
# ├── species1/
# │   ├── config.yaml
# │   ├── sample_metadata.csv          # Required for R analysis
# │   ├── samples/
# │   │   ├── sample1/
# │   │   └── sample2/
# │   └── reference/
# └── species2/
#     ├── config.yaml
#     ├── sample_metadata.csv          # Required for R analysis
#     ├── samples/
#     └── reference/
####################################################################################

set -euo pipefail

##### Version Information #####
PIPELINE_VERSION="2.3.0"
PIPELINE_NAME="IISAGE ATAC-seq Pipeline"

##### Global Configuration from Environment #####

# Base directory containing all species folders (override with IISAGE_BASE_DIR env var)
BASE_DIR="${IISAGE_BASE_DIR:-/mnt/scratch/${USER}/atac_seq_projects}"

# Working directory for SLURM jobs (override with IISAGE_WORK_DIR env var)
WORK_DIR="${IISAGE_WORK_DIR:-/mnt/scratch/${USER}}"

# Email for notifications (override with IISAGE_EMAIL env var)
EMAIL_ADDRESS="${IISAGE_EMAIL:-${USER}@msu.edu}"

# Dry run mode (set DRY_RUN=true to preview without executing)
DRY_RUN="${DRY_RUN:-false}"

# Processing parameters (can be overridden by YAML)
DEFAULT_THREADS=16
DEFAULT_MAX_INSERT_SIZE=2000
DEFAULT_SIZE_SELECTION_CUTOFF=100
DEFAULT_MIN_MAPQ=5
DEFAULT_TRIM_SETTINGS="ILLUMINACLIP:{ADAPTERS_PATH}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:25"
DEFAULT_BOWTIE2_OPTS="--very-sensitive --no-mixed --no-discordant --dovetail"

# DeepTools configurable parameters
DEFAULT_TSS_UPSTREAM=2000
DEFAULT_TSS_DOWNSTREAM=2000
DEFAULT_BIGWIG_BINSIZE=10
DEFAULT_BIGWIG_NORMALIZE="RPKM"

# Disk space requirements (in GB)
MIN_FREE_SPACE_GB=100

# Pipeline control
# Changing these parameters from true/false will change which steps are performed
RUN_TRIMMING=true
RUN_ALIGNMENT=true
RUN_PEAK_CALLING=true
RUN_DEEPTOOLS=true
CLEANUP_INTERMEDIATE=false

# Email notifications
SEND_EMAIL_UPDATES=true

# Resume mode - skip completed stages based on checkpoint files
RESUME_MODE=true

##### R Analysis Configuration #####

# R analysis control
RUN_R_ANALYSIS="${RUN_R_ANALYSIS:-true}"           # Enable/disable R analysis
R_ANALYSIS_REQUIRED="${R_ANALYSIS_REQUIRED:-false}" # Fail if R analysis fails
R_SKIP_MOTIFS="${R_SKIP_MOTIFS:-false}"            # Skip HOMER motif analysis
R_SKIP_ANNOTATION="${R_SKIP_ANNOTATION:-false}"    # Skip peak annotation

# R environment
R_SCRIPT_PATH="${R_SCRIPT_PATH:-$(dirname "$0")/generate_atac_figs_v2.R}"
R_MODULE="${R_MODULE:-R/4.3.0}"                    # R module name for SLURM
R_MIN_SAMPLES=2                                    # Minimum samples for R analysis

# Sample metadata requirements
REQUIRE_SAMPLE_METADATA="${REQUIRE_SAMPLE_METADATA:-true}"  # Require sample_metadata.csv

##### Functions #####

# Create log messages
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log_stage() {
    echo ""
    echo "==========================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] STAGE: $1"
    echo "==========================================="

    if [[ "$SEND_EMAIL_UPDATES" == "true" ]] && [[ "$DRY_RUN" == "false" ]]; then
        send_email "ATAC-seq Pipeline ($CURRENT_SPECIES): $1 Started" \
                   "Stage '$1' has begun for species $CURRENT_SPECIES at $(date '+%Y-%m-%d %H:%M:%S')\nJob ID: $SLURM_JOB_ID\nNode: $SLURM_NODELIST\nVersion: $PIPELINE_VERSION"
    fi
}

# Send emails with better error handling
send_email() {
    if [[ "$SEND_EMAIL_UPDATES" == "true" ]] && [[ "$DRY_RUN" == "false" ]]; then
        if command -v mail >/dev/null 2>&1; then
            local subject="$1"
            local message="$2"
            if ! echo -e "$message" | mail -s "$subject" "$EMAIL_ADDRESS" 2>&1; then
                log_message "WARNING: Failed to send email notification"
            fi
        else
            log_message "WARNING: 'mail' command not found, skipping email notification"
        fi
    fi
}

# Check if files exist
check_file_exists() {
    if [[ ! -f "$1" ]]; then
        log_message "ERROR: Required file not found: $1"
        exit 1
    fi
}

# Check if directories exist
check_dir_exists() {
    if [[ ! -d "$1" ]]; then
        log_message "ERROR: Required directory not found: $1"
        exit 1
    fi
}

# Check available disk space
check_disk_space() {
    local path="$1"
    local required_gb="${2:-$MIN_FREE_SPACE_GB}"

    if [[ "$DRY_RUN" == "true" ]]; then
        log_message "[DRY RUN] Would check disk space on $path (require ${required_gb}GB)"
        return 0
    fi

    local available_kb=$(df -k "$path" | tail -1 | awk '{print $4}')
    local available_gb=$((available_kb / 1024 / 1024))

    log_message "Available disk space: ${available_gb}GB (required: ${required_gb}GB)"

    if [[ $available_gb -lt $required_gb ]]; then
        log_message "ERROR: Insufficient disk space. Available: ${available_gb}GB, Required: ${required_gb}GB"
        exit 1
    fi
}

# Safer YAML parser without eval
parse_yaml_safe() {
    local yaml_file="$1"
    local prefix="$2"

    # Parse YAML and output key=value pairs
    # This version avoids eval for security
    awk -F': ' -v prefix="$prefix" '
    /^[[:space:]]*#/ { next }
    /^[[:space:]]*$/ { next }
    /^[^[:space:]]/ {
        gsub(/[[:space:]]*$/, "", $1)
        if (NF >= 2) {
            # Remove leading/trailing whitespace and quotes
            gsub(/^[[:space:]]*/, "", $2)
            gsub(/[[:space:]]*$/, "", $2)
            gsub(/^["'"'"']|["'"'"']$/, "", $2)
            # Sanitize key name (only allow alphanumeric and underscore)
            key = $1
            gsub(/[^A-Za-z0-9_]/, "_", key)
            print prefix key "=" $2
        }
    }
    ' "$yaml_file"
}

# Set variable safely from key=value string
set_var_safe() {
    local assignment="$1"

    if [[ "$assignment" =~ ^([A-Z_][A-Z0-9_]*)=(.*)$ ]]; then
        local var_name="${BASH_REMATCH[1]}"
        local var_value="${BASH_REMATCH[2]}"

        # Export the variable safely
        export "$var_name"="$var_value"
        log_message "Config: $var_name=$var_value"
    fi
}

# Validate genome size format
validate_genome_size() {
    local size="$1"

    # Check if it's a number or scientific notation (e.g., 2.7e9)
    if [[ "$size" =~ ^[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]+)?$ ]]; then
        return 0
    fi

    log_message "ERROR: Invalid genome size format: $size"
    return 1
}

# Load species configuration from YAML file
load_species_config() {
    local species_dir="$1"
    local config_file="$species_dir/config.yaml"

    if [[ ! -f "$config_file" ]]; then
        log_message "ERROR: Config file not found: $config_file"
        return 1
    fi

    log_message "Loading configuration from: $config_file"

    # Parse YAML and set variables safely (no eval)
    while IFS= read -r line; do
        if [[ "$line" =~ ^[A-Z_]+=.* ]]; then
            set_var_safe "$line"
        fi
    done < <(parse_yaml_safe "$config_file" "")

    # Set species-specific paths
    export SPECIES_DIR="$species_dir"
    export SAMPLES_DIR="$species_dir/samples"
    export PROCESSED_DIR="$species_dir/processed_files"
    export PEAKS_DIR="$species_dir/peaks_output"
    export DEEPTOOLS_DIR="$species_dir/deeptools_output"
    export CHECKPOINT_DIR="$species_dir/.checkpoints"

    # Validate required config variables
    local required_vars=("GENOME_SIZE" "REF_GENOME_FASTA" "SPECIES_NAME")
    for var in "${required_vars[@]}"; do
        if [[ -z "${!var:-}" ]]; then
            log_message "ERROR: Required variable $var not found in config.yaml"
            return 1
        fi
    done

    # Validate genome size format
    if ! validate_genome_size "$GENOME_SIZE"; then
        return 1
    fi

    # Set optional variables with defaults
    export THREADS="${THREADS:-$DEFAULT_THREADS}"
    export MAX_INSERT_SIZE="${MAX_INSERT_SIZE:-$DEFAULT_MAX_INSERT_SIZE}"
    export SIZE_SELECTION_CUTOFF="${SIZE_SELECTION_CUTOFF:-$DEFAULT_SIZE_SELECTION_CUTOFF}"
    export MIN_MAPQ="${MIN_MAPQ:-$DEFAULT_MIN_MAPQ}"
    export BOWTIE2_OPTS="${BOWTIE2_OPTS:-$DEFAULT_BOWTIE2_OPTS}"
    export TSS_UPSTREAM="${TSS_UPSTREAM:-$DEFAULT_TSS_UPSTREAM}"
    export TSS_DOWNSTREAM="${TSS_DOWNSTREAM:-$DEFAULT_TSS_DOWNSTREAM}"
    export BIGWIG_BINSIZE="${BIGWIG_BINSIZE:-$DEFAULT_BIGWIG_BINSIZE}"
    export BIGWIG_NORMALIZE="${BIGWIG_NORMALIZE:-$DEFAULT_BIGWIG_NORMALIZE}"

    # Handle relative paths in config
    if [[ ! "$REF_GENOME_FASTA" = /* ]]; then
        export REF_GENOME_FASTA="$species_dir/$REF_GENOME_FASTA"
    fi
    if [[ ! "${ADAPTERS_PATH:-}" = /* ]] && [[ -n "${ADAPTERS_PATH:-}" ]]; then
        export ADAPTERS_PATH="$species_dir/$ADAPTERS_PATH"
    fi
    if [[ ! "${REF_GENOME:-}" = /* ]] && [[ -n "${REF_GENOME:-}" ]]; then
        export REF_GENOME="$species_dir/$REF_GENOME"
    else
        # Default genome index location
        export REF_GENOME="$species_dir/reference/genome_index"
    fi

    # Set TRIM_SETTINGS with actual adapters path
    if [[ -n "${ADAPTERS_PATH:-}" ]]; then
        export TRIM_SETTINGS="${TRIM_SETTINGS:-$DEFAULT_TRIM_SETTINGS}"
        export TRIM_SETTINGS="${TRIM_SETTINGS//\{ADAPTERS_PATH\}/$ADAPTERS_PATH}"
    fi

    return 0
}

# Load modules for Stage 1: Trimming & Alignment
load_trimming_alignment_modules() {
    log_message "Loading Stage 1 modules (Trimming & Alignment)"

    if [[ "$DRY_RUN" == "true" ]]; then
        log_message "[DRY RUN] Stage 1 would load: Trimmomatic, Bowtie2/2.5.1-GCC-12.2.0, SAMtools/1.17-GCC-12.2.0"
        return 0
    fi

    # Clean environment completely
    module purge

    # Load Stage 1 modules (compatible with GCC 12.2.0 / foss-2022b)
    local modules=(
        "Trimmomatic"
        "Bowtie2/2.5.1-GCC-12.2.0"
        "SAMtools/1.17-GCC-12.2.0"
    )

    local failed_modules=()

    for mod in "${modules[@]}"; do
        if module load "$mod" 2>&1; then
            log_message "✓ Loaded: $mod"
        else
            log_message "✗ Failed to load module: $mod"
            failed_modules+=("$mod")
        fi
    done

    if [[ ${#failed_modules[@]} -gt 0 ]]; then
        log_message "ERROR: Failed to load Stage 1 modules: ${failed_modules[*]}"
        log_message "Available modules can be checked with: module spider <module_name>"
        exit 1
    fi

    # Verify key commands/variables exist
    # Trimmomatic is run via java -jar, so check for environment variable
    if [[ -z "${EBROOTTRIMMOMATIC:-}" ]]; then
        log_message "ERROR: EBROOTTRIMMOMATIC not set after loading Trimmomatic module"
        exit 1
    fi

    # Check for standard commands
    local commands=("bowtie2" "samtools")
    for cmd in "${commands[@]}"; do
        if ! command -v "$cmd" >/dev/null 2>&1; then
            log_message "ERROR: Command '$cmd' not found after loading modules"
            exit 1
        fi
    done

    log_message "✓ Stage 1 modules loaded successfully"
}

# Load modules for Stage 2: Peak Calling
load_peakcalling_modules() {
    log_message "Loading Stage 2 modules (Peak Calling)"

    if [[ "$DRY_RUN" == "true" ]]; then
        log_message "[DRY RUN] Stage 2 would load: SAMtools/1.17-GCC-12.2.0, MACS2/2.2.9.1-foss-2022b"
        return 0
    fi

    # Clean environment completely
    module purge

    # Load Stage 2 modules (MACS2 requires foss-2022b)
    local modules=(
        "SAMtools/1.17-GCC-12.2.0"
        "MACS2/2.2.9.1-foss-2022b"
    )

    local failed_modules=()

    for mod in "${modules[@]}"; do
        if module load "$mod" 2>&1; then
            log_message "✓ Loaded: $mod"
        else
            log_message "✗ Failed to load module: $mod"
            failed_modules+=("$mod")
        fi
    done

    if [[ ${#failed_modules[@]} -gt 0 ]]; then
        log_message "ERROR: Failed to load Stage 2 modules: ${failed_modules[*]}"
        log_message "Available modules can be checked with: module spider <module_name>"
        exit 1
    fi

    # Verify key commands exist
    local commands=("samtools" "macs2")
    for cmd in "${commands[@]}"; do
        if ! command -v "$cmd" >/dev/null 2>&1; then
            log_message "ERROR: Command '$cmd' not found after loading modules"
            exit 1
        fi
    done

    log_message "✓ Stage 2 modules loaded successfully"
}

# Load modules for Stage 3: DeepTools Visualization
load_deeptools_modules() {
    log_message "Loading Stage 3 modules (DeepTools Visualization)"

    if [[ "$DRY_RUN" == "true" ]]; then
        log_message "[DRY RUN] Stage 3 would load: SAMtools, deepTools/3.5.5-foss-2023a"
        return 0
    fi

    # Clean environment completely
    module purge

    # Load Stage 3 modules (deepTools requires foss-2023a)
    local modules=(
        "SAMtools"  # Will load version compatible with foss-2023a
        "deepTools/3.5.5-foss-2023a"
    )

    local failed_modules=()

    for mod in "${modules[@]}"; do
        if module load "$mod" 2>&1; then
            log_message "✓ Loaded: $mod"
        else
            log_message "✗ Failed to load module: $mod"
            failed_modules+=("$mod")
        fi
    done

    if [[ ${#failed_modules[@]} -gt 0 ]]; then
        log_message "ERROR: Failed to load Stage 3 modules: ${failed_modules[*]}"
        log_message "Available modules can be checked with: module spider <module_name>"
        exit 1
    fi

    # Verify key commands exist
    local commands=("samtools" "bamCoverage" "computeMatrix" "plotHeatmap" "multiBamSummary" "plotCorrelation")
    for cmd in "${commands[@]}"; do
        if ! command -v "$cmd" >/dev/null 2>&1; then
            log_message "ERROR: Command '$cmd' not found after loading modules"
            exit 1
        fi
    done

    log_message "✓ Stage 3 modules loaded successfully"
}

# Record software versions for reproducibility
record_versions() {
    local version_file="$1"

    if [[ "$DRY_RUN" == "true" ]]; then
        log_message "[DRY RUN] Would record software versions to $version_file"
        return 0
    fi

    {
        echo "=== Software Versions ==="
        echo "Pipeline: $PIPELINE_NAME v$PIPELINE_VERSION"
        echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "Host: $(hostname)"
        echo "User: $USER"
        echo ""

        echo "--- Core Tools ---"
        echo "Note: Versions recorded from Stage 1 modules"
        echo "Trimmomatic: See module list"
        echo "Bowtie2: See module list"
        echo "SAMtools: See module list"
        echo "MACS2: See module list (loaded in Stage 2)"
        echo "deepTools: See module list (loaded in Stage 3)"

        echo ""
        echo "--- Environment ---"
        echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-N/A}"
        echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-N/A}"
        echo "CPUs: $THREADS"
        echo "Genome: $REF_GENOME_FASTA"
        echo "Genome Size: $GENOME_SIZE"

    } > "$version_file"

    log_message "Software versions recorded to: $version_file"
}

# Checkpoint functions for resume capability
create_checkpoint() {
    local checkpoint_name="$1"
    local checkpoint_file="$CHECKPOINT_DIR/${checkpoint_name}.done"

    if [[ "$DRY_RUN" == "true" ]]; then
        log_message "[DRY RUN] Would create checkpoint: $checkpoint_name"
        return 0
    fi

    mkdir -p "$CHECKPOINT_DIR"
    echo "$(date '+%Y-%m-%d %H:%M:%S')" > "$checkpoint_file"
    log_message "Checkpoint created: $checkpoint_name"
}

check_checkpoint() {
    local checkpoint_name="$1"
    local checkpoint_file="$CHECKPOINT_DIR/${checkpoint_name}.done"

    if [[ "$RESUME_MODE" == "true" ]] && [[ -f "$checkpoint_file" ]]; then
        log_message "Checkpoint exists: $checkpoint_name (created: $(cat "$checkpoint_file"))"
        return 0
    else
        return 1
    fi
}

# Deploy species job (submitting an individual SLURM job for each species)
deploy_species_job() {
    local species_name="$1"
    local species_dir="$BASE_DIR/$species_name"

    if [[ ! -d "$species_dir" ]]; then
        log_message "WARNING: Species directory not found: $species_dir"
        return 1
    fi

    log_message "Deploying job for species: $species_name"

    if [[ "$DRY_RUN" == "true" ]]; then
        log_message "[DRY RUN] Would submit SLURM job for $species_name"
        return 0
    fi

    # Count samples for array job
    local sample_count=$(find "$species_dir/samples" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)

    if [[ $sample_count -gt 0 ]]; then
        # Submit job for this species
        sbatch --export=SPECIES_NAME="$species_name",IISAGE_BASE_DIR="$BASE_DIR",IISAGE_EMAIL="$EMAIL_ADDRESS" \
               --job-name="atac-${species_name}" \
               --output="atac-${species_name}-%j.SLURMout" \
               --error="atac-${species_name}-%j.SLURMerr" \
               "$0"

        log_message "Job submitted for $species_name ($sample_count samples)"
    else
        log_message "WARNING: No samples found for $species_name, skipping"
    fi
}

# Remove intermediate files
cleanup_files() {
    local files=("$@")
    if [[ "$CLEANUP_INTERMEDIATE" == "true" ]]; then
        for file in "${files[@]}"; do
            if [[ -f "$file" ]]; then
                if [[ "$DRY_RUN" == "true" ]]; then
                    log_message "[DRY RUN] Would clean up: $(basename "$file")"
                else
                    log_message "Cleaning up: $(basename "$file")"
                    rm "$file"
                fi
            fi
        done
    fi
}

# Execute command with error checking
execute_cmd() {
    local cmd="$1"
    local error_msg="${2:-Command failed}"

    if [[ "$DRY_RUN" == "true" ]]; then
        log_message "[DRY RUN] Would execute: $cmd"
        return 0
    fi

    log_message "Executing: $cmd"
    if ! eval "$cmd"; then
        log_message "ERROR: $error_msg"
        log_message "Failed command: $cmd"
        exit 1
    fi
}

# Load R module and validate R environment
load_r_environment() {
    if [[ "$RUN_R_ANALYSIS" != "true" ]]; then
        return 0
    fi

    log_message "Validating R environment for analysis"

    if [[ "$DRY_RUN" == "true" ]]; then
        log_message "[DRY RUN] Would load R module: $R_MODULE"
        log_message "[DRY RUN] Would validate R script: $R_SCRIPT_PATH"
        return 0
    fi

    # Clean environment and load R module
    module purge

    if ! module load "$R_MODULE" 2>&1; then
        log_message "WARNING: Failed to load R module: $R_MODULE"
        if [[ "$R_ANALYSIS_REQUIRED" == "true" ]]; then
            log_message "ERROR: R analysis required but R module unavailable"
            return 1
        else
            log_message "R analysis will be skipped"
            export RUN_R_ANALYSIS=false
            return 0
        fi
    fi

    # Verify R is available
    if ! command -v Rscript >/dev/null 2>&1; then
        log_message "WARNING: Rscript not found after loading R module"
        export RUN_R_ANALYSIS=false
        return 0
    fi

    # Check if R script exists
    if [[ ! -f "$R_SCRIPT_PATH" ]]; then
        log_message "WARNING: R script not found: $R_SCRIPT_PATH"
        if [[ "$R_ANALYSIS_REQUIRED" == "true" ]]; then
            log_message "ERROR: R analysis required but script not found"
            return 1
        else
            export RUN_R_ANALYSIS=false
            return 0
        fi
    fi

    log_message "R environment validated successfully"
    log_message "  R version: $(Rscript --version 2>&1 | head -1)"
    log_message "  R script: $R_SCRIPT_PATH"
    return 0
}

# Validate sample metadata file
validate_sample_metadata() {
    local species_dir="$1"
    local metadata_file="$species_dir/sample_metadata.csv"

    log_message "Validating sample metadata file"

    # Check file exists
    if [[ ! -f "$metadata_file" ]]; then
        log_message "ERROR: sample_metadata.csv not found: $metadata_file"
        log_message "Create a file with columns: sampID,condition,replicate"
        log_message "Example:"
        log_message "  sampID,condition,replicate"
        log_message "  sample1,control,1"
        log_message "  sample2,treatment,1"
        return 1
    fi

    # Check file is not empty
    if [[ ! -s "$metadata_file" ]]; then
        log_message "ERROR: Metadata file is empty: $metadata_file"
        return 1
    fi

    # Check header
    local header=$(head -n1 "$metadata_file")
    if [[ "$header" != "sampID,condition,replicate" ]]; then
        log_message "WARNING: Unexpected metadata header: $header"
        log_message "Expected: sampID,condition,replicate"
    fi

    # Count samples
    local sample_count=$(tail -n +2 "$metadata_file" | grep -v '^[[:space:]]*$' | wc -l)
    log_message "Metadata contains $sample_count samples"

    if [[ $sample_count -lt $R_MIN_SAMPLES ]]; then
        log_message "WARNING: Only $sample_count samples (minimum $R_MIN_SAMPLES recommended)"
    fi

    if [[ $sample_count -eq 0 ]]; then
        log_message "ERROR: No samples found in metadata file"
        return 1
    fi

    log_message "Sample metadata validation successful"
    return 0
}

# Run R analysis on completed peak data
run_r_analysis() {
    local species_dir="$1"
    local base_dir="$(dirname "$species_dir")"

    log_message "Running R analysis"

    if [[ "$DRY_RUN" == "true" ]]; then
        local r_cmd="Rscript '$R_SCRIPT_PATH' --base_dir='$base_dir'"
        if [[ "$R_SKIP_MOTIFS" == "true" ]]; then
            r_cmd="$r_cmd --no-motifs"
        fi
        if [[ "$R_SKIP_ANNOTATION" == "true" ]]; then
            r_cmd="$r_cmd --no-annotation"
        fi
        log_message "[DRY RUN] Would execute: $r_cmd"
        return 0
    fi

    # Build R command
    local r_cmd="Rscript '$R_SCRIPT_PATH' --base_dir='$base_dir'"

    if [[ "$R_SKIP_MOTIFS" == "true" ]]; then
        r_cmd="$r_cmd --no-motifs"
    fi

    if [[ "$R_SKIP_ANNOTATION" == "true" ]]; then
        r_cmd="$r_cmd --no-annotation"
    fi

    log_message "Executing: $r_cmd"

    # Run R analysis with error handling
    local r_output_log="$PROCESSED_DIR/r_analysis_output.log"
    local r_error_log="$PROCESSED_DIR/r_analysis_error.log"

    if eval "$r_cmd" > "$r_output_log" 2> "$r_error_log"; then
        log_message "R analysis completed successfully"
        log_message "Output logged to: $r_output_log"

        # Check for generated figures
        local figures_dir="$species_dir/figures"
        if [[ -d "$figures_dir" ]]; then
            local figure_count=$(find "$figures_dir" -type f \( -name "*.png" -o -name "*.pdf" \) 2>/dev/null | wc -l)
            local csv_count=$(find "$figures_dir" -type f -name "*.csv" 2>/dev/null | wc -l)
            log_message "Generated $figure_count figures and $csv_count tables in: $figures_dir"
        fi

        return 0
    else
        log_message "WARNING: R analysis failed"
        log_message "Error log: $r_error_log"

        # Show last 20 lines of error
        if [[ -f "$r_error_log" ]]; then
            log_message "Last 20 lines of error log:"
            tail -n 20 "$r_error_log" | while IFS= read -r line; do
                log_message "  $line"
            done
        fi

        if [[ "$R_ANALYSIS_REQUIRED" == "true" ]]; then
            log_message "ERROR: R analysis required but failed"
            return 1
        else
            log_message "Continuing despite R analysis failure"
            return 0
        fi
    fi
}

# Build or validate Bowtie2 index with smart array job handling
setup_bowtie2_index() {
    local ref_fasta="$1"
    local ref_genome="$2"
    local lock_dir="$CHECKPOINT_DIR/bowtie2_index.lock"
    local max_wait_seconds=7200  # 2 hours
    local poll_interval=30        # 30 seconds
    local stale_lock_hours=3      # Consider lock stale after 3 hours

    # Check if all index files exist
    local index_complete=true
    for ext in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
        if [[ ! -f "${ref_genome}.${ext}" ]]; then
            index_complete=false
            break
        fi
    done

    if [[ "$index_complete" == "true" ]]; then
        log_message "Bowtie2 index found and complete for $SPECIES_NAME"
        return 0
    fi

    if check_checkpoint "bowtie2_index"; then
        log_message "Bowtie2 index checkpoint exists, skipping rebuild"
        return 0
    fi

    # Array job handling: Use lock mechanism to coordinate index building
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log_message "Running in array mode (task ${SLURM_ARRAY_TASK_ID})"

        # Check for stale lock (another job crashed)
        if [[ -d "$lock_dir" ]]; then
            local lock_age_seconds=$(( $(date +%s) - $(stat -c %Y "$lock_dir" 2>/dev/null || echo 0) ))
            local stale_threshold=$((stale_lock_hours * 3600))

            if [[ $lock_age_seconds -gt $stale_threshold ]]; then
                log_message "WARNING: Removing stale lock (${lock_age_seconds}s old, threshold: ${stale_threshold}s)"
                rmdir "$lock_dir" 2>/dev/null || true
            fi
        fi

        # Try to acquire lock (atomic operation)
        if mkdir "$lock_dir" 2>/dev/null; then
            # This task acquired the lock - it will build the index
            log_message "Array task ${SLURM_ARRAY_TASK_ID}: Acquired lock, building index for all tasks"

            trap "rmdir '$lock_dir' 2>/dev/null || true" EXIT ERR

            # Build index
            log_message "Building Bowtie2 index for $SPECIES_NAME..."
            mkdir -p "$(dirname "$ref_genome")"

            execute_cmd "bowtie2-build '$ref_fasta' '$ref_genome'" \
                        "Failed to build Bowtie2 index"

            # Validate index was created
            if [[ ! -f "${ref_genome}.1.bt2" ]]; then
                log_message "ERROR: Bowtie2 index build completed but index files not found"
                rmdir "$lock_dir" 2>/dev/null || true
                exit 1
            fi

            create_checkpoint "bowtie2_index"
            log_message "Index created successfully by array task ${SLURM_ARRAY_TASK_ID}"

            # Release lock
            rmdir "$lock_dir" 2>/dev/null || true
            trap - EXIT ERR

        else
            # Lock exists - another task is building the index
            log_message "Array task ${SLURM_ARRAY_TASK_ID}: Another task is building index, waiting..."

            local waited=0
            local max_attempts=$((max_wait_seconds / poll_interval))
            local attempt=0

            while [[ $waited -lt $max_wait_seconds ]]; do
                ((attempt++))

                # Check if index is complete
                local index_ready=true
                for ext in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
                    if [[ ! -f "${ref_genome}.${ext}" ]]; then
                        index_ready=false
                        break
                    fi
                done

                if [[ "$index_ready" == "true" ]] && check_checkpoint "bowtie2_index"; then
                    log_message "Array task ${SLURM_ARRAY_TASK_ID}: Index ready after ${waited}s, proceeding"
                    return 0
                fi

                # Log progress every 5 minutes
                if (( attempt % 10 == 0 )); then
                    log_message "Array task ${SLURM_ARRAY_TASK_ID}: Still waiting for index (${waited}s / ${max_wait_seconds}s)"
                fi

                sleep $poll_interval
                waited=$((waited + poll_interval))
            done

            # Timeout reached
            log_message "ERROR: Timeout waiting for index build after ${max_wait_seconds}s"
            log_message "  Possible causes:"
            log_message "    - Index build taking longer than expected (large genome)"
            log_message "    - Builder task failed or was cancelled"
            log_message "  Solutions:"
            log_message "    - Check other array task logs for build errors"
            log_message "    - Build index separately first (non-array job)"
            log_message "    - Increase timeout (MAX_WAIT_SECONDS in script)"
            exit 1
        fi
    else
        # Not an array job - build normally
        log_message "Building Bowtie2 index for $SPECIES_NAME..."
        mkdir -p "$(dirname "$ref_genome")"

        execute_cmd "bowtie2-build '$ref_fasta' '$ref_genome'" \
                    "Failed to build Bowtie2 index"

        # Validate index was created
        if [[ ! -f "${ref_genome}.1.bt2" ]]; then
            log_message "ERROR: Bowtie2 index build completed but index files not found"
            exit 1
        fi

        create_checkpoint "bowtie2_index"
        log_message "Index created successfully"
    fi
}

# Get list of sample directories
get_sample_list() {
    local samples_dir="$1"
    find "$samples_dir" -mindepth 1 -maxdepth 1 -type d | sort
}

# Process single sample (for array jobs)
process_sample_array() {
    local sample_dir="$1"
    local sample_name=$(basename "$sample_dir")

    log_message "Processing sample (array mode): $sample_name"

    # Check for input files
    local r1_files=("$sample_dir"/*_1.fq.gz)
    local r2_files=("$sample_dir"/*_2.fq.gz)

    if [[ ! -e "${r1_files[0]}" ]] || [[ ! -e "${r2_files[0]}" ]]; then
        log_message "WARNING: No paired FASTQ files found in $sample_dir, skipping"
        return 1
    fi

    # Define file paths
    local bam_file="$PROCESSED_DIR/${sample_name}_sorted.bam"

    # Check if already processed
    if check_checkpoint "sample_${sample_name}_aligned"; then
        log_message "Sample $sample_name already processed, skipping"
        return 0
    fi

    if [[ "$RUN_TRIMMING" == "true" ]] && [[ "$RUN_ALIGNMENT" == "true" ]]; then
        # File paths
        local merged_r1="$PROCESSED_DIR/${sample_name}_R1_merged.fastq"
        local merged_r2="$PROCESSED_DIR/${sample_name}_R2_merged.fastq"
        local paired_r1="$PROCESSED_DIR/${sample_name}_R1_paired.fastq"
        local paired_r2="$PROCESSED_DIR/${sample_name}_R2_paired.fastq"
        local unpaired_r1="$PROCESSED_DIR/${sample_name}_R1_unpaired.fastq"
        local unpaired_r2="$PROCESSED_DIR/${sample_name}_R2_unpaired.fastq"
        local sam_file="$PROCESSED_DIR/${sample_name}_aligned.sam"

        # Merge FASTQ files
        log_message "Merging FASTQ files for $sample_name"
        execute_cmd "zcat '$sample_dir'/*_1.fq.gz > '$merged_r1'" \
                    "Failed to merge R1 FASTQ files"
        execute_cmd "zcat '$sample_dir'/*_2.fq.gz > '$merged_r2'" \
                    "Failed to merge R2 FASTQ files"

        # Trim adapters
        log_message "Trimming adapters for $sample_name"
        execute_cmd "java -jar \"\$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar\" PE -threads $THREADS \
            '$merged_r1' '$merged_r2' \
            '$paired_r1' '$unpaired_r1' '$paired_r2' '$unpaired_r2' \
            $TRIM_SETTINGS" \
            "Failed to trim adapters for $sample_name"

        cleanup_files "$merged_r1" "$merged_r2"

        # Alignment
        log_message "Aligning reads for $sample_name"
        execute_cmd "bowtie2 -x '$REF_GENOME' -1 '$paired_r1' -2 '$paired_r2' \
                -S '$sam_file' $BOWTIE2_OPTS \
                -X $MAX_INSERT_SIZE -k 2 -p $THREADS" \
                "Failed to align reads for $sample_name"

        # Convert to BAM, sort and index
        log_message "Converting to BAM and sorting for $sample_name"
        execute_cmd "samtools view -bS '$sam_file' | samtools sort -@ $THREADS -o '$bam_file'" \
                    "Failed to convert/sort BAM for $sample_name"

        execute_cmd "samtools index '$bam_file'" \
                    "Failed to index BAM for $sample_name"

        # Validate BAM file
        if [[ "$DRY_RUN" == "false" ]]; then
            if ! samtools quickcheck "$bam_file" 2>/dev/null; then
                log_message "ERROR: BAM file validation failed for $sample_name"
                exit 1
            fi
        fi

        cleanup_files "$paired_r1" "$paired_r2" "$unpaired_r1" "$unpaired_r2" "$sam_file"

        create_checkpoint "sample_${sample_name}_aligned"
        log_message "Successfully processed sample: $sample_name"
    fi
}

##### Main Pipeline Logic #####

# Print header
log_message "========================================="
log_message "$PIPELINE_NAME v$PIPELINE_VERSION"
log_message "========================================="

if [[ "$DRY_RUN" == "true" ]]; then
    log_message "*** DRY RUN MODE - No commands will be executed ***"
fi

# Change to working directory
if [[ "$DRY_RUN" == "false" ]]; then
    cd "$WORK_DIR"
fi

# Check if this is a deployment run
if [[ "${DEPLOY_ALL:-false}" == "true" ]]; then
    log_stage "Multi-Species Pipeline Deployment"

    check_dir_exists "$BASE_DIR"

    # Find all species directories (those containing config.yaml)
    species_count=0
    for species_dir in "$BASE_DIR"/*/; do
        if [[ -f "$species_dir/config.yaml" ]]; then
            species_name=$(basename "$species_dir")
            deploy_species_job "$species_name"
            ((species_count++))
        fi
    done

    log_message "Deployed $species_count species jobs"
    send_email "ATAC-seq Multi-Species Deployment Complete" \
               "Deployed $species_count species jobs at $(date '+%Y-%m-%d %H:%M:%S')\nPipeline Version: $PIPELINE_VERSION"
    exit 0
fi

# Single species processing
if [[ -z "${SPECIES_NAME:-}" ]]; then
    log_message "ERROR: SPECIES_NAME not provided. Use --export=SPECIES_NAME=species_name"
    exit 1
fi

export CURRENT_SPECIES="$SPECIES_NAME"
SPECIES_DIR="$BASE_DIR/$SPECIES_NAME"

log_stage "Single Species Pipeline: $SPECIES_NAME"

# Load species configuration
if ! load_species_config "$SPECIES_DIR"; then
    log_message "ERROR: Failed to load configuration for $SPECIES_NAME"
    exit 1
fi

# Check disk space before proceeding
check_disk_space "$SPECIES_DIR"

# Load R environment if needed
if ! load_r_environment; then
    if [[ "$R_ANALYSIS_REQUIRED" == "true" ]]; then
        log_message "ERROR: R environment required but unavailable"
        exit 1
    fi
fi

# Validate paths
check_dir_exists "$SAMPLES_DIR"
check_file_exists "$REF_GENOME_FASTA"

# Create output directories
if [[ "$DRY_RUN" == "false" ]]; then
    mkdir -p "$PROCESSED_DIR" "$PEAKS_DIR" "$DEEPTOOLS_DIR" "$CHECKPOINT_DIR"
fi

# Record software versions
VERSION_FILE="$PROCESSED_DIR/software_versions.txt"
record_versions "$VERSION_FILE"

# Log configuration summary
log_message "Configuration summary for $SPECIES_NAME:"
log_message "  Genome size: $GENOME_SIZE"
log_message "  Reference: $REF_GENOME_FASTA"
log_message "  Samples directory: $SAMPLES_DIR"
log_message "  Threads: $THREADS"
log_message "  TSS region: -${TSS_UPSTREAM}bp to +${TSS_DOWNSTREAM}bp"
log_message "  BigWig bin size: ${BIGWIG_BINSIZE}bp"
log_message "  BigWig normalization: $BIGWIG_NORMALIZE"
log_message "  Resume mode: $RESUME_MODE"

##### STAGE 1: Trimming & Alignment #####

if [[ "$RUN_TRIMMING" == "true" ]] || [[ "$RUN_ALIGNMENT" == "true" ]]; then
    log_stage "Quality Control, Trimming & Alignment"

    # Load Stage 1 modules
    load_trimming_alignment_modules

    # Setup Bowtie2 index
    setup_bowtie2_index "$REF_GENOME_FASTA" "$REF_GENOME"

    # Get list of samples
    mapfile -t sample_dirs < <(get_sample_list "$SAMPLES_DIR")
    total_samples=${#sample_dirs[@]}

    log_message "Found $total_samples samples for $SPECIES_NAME"

    # Check if running as SLURM array job
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        # Array mode - process single sample
        array_index=$((SLURM_ARRAY_TASK_ID - 1))
        if [[ $array_index -lt $total_samples ]]; then
            sample_dir="${sample_dirs[$array_index]}"
            process_sample_array "$sample_dir"
        else
            log_message "ERROR: Array task ID $SLURM_ARRAY_TASK_ID exceeds sample count $total_samples"
            exit 1
        fi
    else
        # Sequential mode - process all samples
        current_sample=0
        for sample_dir in "${sample_dirs[@]}"; do
            ((current_sample++))
            sample_name=$(basename "$sample_dir")

            log_message "Processing sample $current_sample/$total_samples: $sample_name ($SPECIES_NAME)"
            process_sample_array "$sample_dir"
        done
    fi
fi

# If running in array mode, exit here (peak calling and visualization happen in main job)
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log_message "Array job completed for task ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

##### STAGE 2: Peak Calling with Size Selection #####

if [[ "$RUN_PEAK_CALLING" == "true" ]]; then
    log_stage "Peak Calling with Size Selection"

    # Load Stage 2 modules
    load_peakcalling_modules

    for bam_file in "$PROCESSED_DIR"/*_sorted.bam; do
        if [[ ! -f "$bam_file" ]]; then
            log_message "No BAM files found in $PROCESSED_DIR"
            break
        fi

        sample_name=$(basename "$bam_file" _sorted.bam)

        # Check checkpoint
        if check_checkpoint "sample_${sample_name}_peaks"; then
            log_message "Peaks already called for $sample_name, skipping"
            continue
        fi

        log_message "Processing peaks for $sample_name ($SPECIES_NAME)"

        # Size selection and filtering
        log_message "Performing size selection and filtering for $sample_name"
        size_selected_bam="${PEAKS_DIR}/${sample_name}_size_selected_${SIZE_SELECTION_CUTOFF}.bam"

        execute_cmd "samtools view -@ $THREADS -h '$bam_file' | \
        awk -v LEN=$SIZE_SELECTION_CUTOFF '{
            if (\$9 <= LEN && \$9 >= -LEN && \$9 != 0 || \$1 ~ /^@/) print \$0
        }' | \
        grep -v '_random\|chrUn\|chrM' | \
        samtools view -@ $THREADS -uh -q $MIN_MAPQ - | \
        samtools sort -@ $THREADS -o '$size_selected_bam'" \
        "Failed size selection for $sample_name"

        execute_cmd "samtools index -@ $THREADS '$size_selected_bam'" \
                    "Failed to index size-selected BAM for $sample_name"

        # Call peaks with MACS2
        log_message "Calling peaks for $sample_name"
        execute_cmd "macs2 callpeak -t '$size_selected_bam' \
                       -f BAM \
                       -g $GENOME_SIZE \
                       --nomodel \
                       --nolambda \
                       --keep-dup all \
                       --call-summits \
                       -n '$sample_name' \
                       --outdir '$PEAKS_DIR'" \
                    "Failed to call peaks for $sample_name"

        create_checkpoint "sample_${sample_name}_peaks"
        log_message "Peaks called for $sample_name"
    done
fi

##### STAGE 3: DeepTools Visualization #####

if [[ "$RUN_DEEPTOOLS" == "true" ]]; then
    log_stage "DeepTools Visualization"

    # Load Stage 3 modules
    load_deeptools_modules

    if check_checkpoint "deeptools_complete"; then
        log_message "DeepTools analysis already completed, skipping"
    else
        # Collect all size-selected BAM files
        size_selected_bams=("$PEAKS_DIR"/*_size_selected_${SIZE_SELECTION_CUTOFF}.bam)

        if [[ ${#size_selected_bams[@]} -eq 0 ]] || [[ ! -f "${size_selected_bams[0]}" ]]; then
            log_message "No size-selected BAM files found for DeepTools analysis"
        else
            log_message "Running DeepTools analysis on ${#size_selected_bams[@]} samples for $SPECIES_NAME"

            # Create sample labels
            sample_labels=""
            for bam in "${size_selected_bams[@]}"; do
                sample_name=$(basename "$bam" "_size_selected_${SIZE_SELECTION_CUTOFF}.bam")
                sample_labels="$sample_labels $sample_name"
            done

            # Look for TSS file in species directory
            tss_bed="$SPECIES_DIR/reference/tss_regions.bed"

            if [[ -f "$tss_bed" ]]; then
                log_message "Computing matrix around TSS regions (±${TSS_UPSTREAM}/${TSS_DOWNSTREAM}bp)"
                execute_cmd "computeMatrix reference-point \
                    --referencePoint TSS \
                    -b $TSS_UPSTREAM -a $TSS_DOWNSTREAM \
                    -R '$tss_bed' \
                    -S ${size_selected_bams[*]} \
                    --skipZeros \
                    -o '$DEEPTOOLS_DIR/matrix_TSS.gz' \
                    --samplesLabel $sample_labels \
                    -p $THREADS" \
                    "Failed to compute TSS matrix"

                # Create heatmap
                log_message "Creating TSS heatmap"
                execute_cmd "plotHeatmap -m '$DEEPTOOLS_DIR/matrix_TSS.gz' \
                           -out '$DEEPTOOLS_DIR/TSS_heatmap.png' \
                           --colorMap viridis \
                           --whatToShow 'heatmap and colorbar'" \
                           "Failed to create TSS heatmap"
            else
                log_message "TSS BED file not found at $tss_bed, skipping TSS analysis"
            fi

            # Compute coverage correlation
            log_message "Computing sample correlation"
            execute_cmd "multiBamSummary bins \
                --bamfiles ${size_selected_bams[*]} \
                --labels $sample_labels \
                -out '$DEEPTOOLS_DIR/multiBamSummary.npz' \
                -p $THREADS" \
                "Failed to compute BAM summary"

            # Plot correlation
            execute_cmd "plotCorrelation \
                -in '$DEEPTOOLS_DIR/multiBamSummary.npz' \
                --corMethod pearson \
                --skipZeros \
                --plotTitle 'Sample Correlation - $SPECIES_NAME' \
                --whatToPlot heatmap \
                --colorMap RdYlBu \
                --plotNumbers \
                -o '$DEEPTOOLS_DIR/correlation_heatmap.png'" \
                "Failed to create correlation plot"

            # Create bigWig files
            log_message "Creating bigWig files for genome browser visualization (${BIGWIG_NORMALIZE} normalized, ${BIGWIG_BINSIZE}bp bins)"
            for bam in "${size_selected_bams[@]}"; do
                sample_name=$(basename "$bam" "_size_selected_${SIZE_SELECTION_CUTOFF}.bam")
                log_message "Creating bigWig for $sample_name"
                execute_cmd "bamCoverage -b '$bam' \
                           -o '$DEEPTOOLS_DIR/${sample_name}.bw' \
                           --normalizeUsing $BIGWIG_NORMALIZE \
                           --binSize $BIGWIG_BINSIZE \
                           -p $THREADS" \
                           "Failed to create bigWig for $sample_name"
            done

            create_checkpoint "deeptools_complete"
        fi
    fi
fi

##### STAGE 4: R Analysis (Downstream Analysis & Visualization) #####

if [[ "$RUN_R_ANALYSIS" == "true" ]]; then
    log_stage "R-based Downstream Analysis"

    if check_checkpoint "r_analysis_complete"; then
        log_message "R analysis already completed, skipping"
    else
        # Check we have samples processed
        peak_count=$(find "$PEAKS_DIR" -name "*_peaks.narrowPeak" -type f 2>/dev/null | wc -l)
        bam_count=$(find "$PEAKS_DIR" -name "*_size_selected_${SIZE_SELECTION_CUTOFF}.bam" -type f 2>/dev/null | wc -l)

        if [[ $peak_count -eq 0 ]] || [[ $bam_count -eq 0 ]]; then
            log_message "WARNING: No peaks or BAM files found for R analysis"
            log_message "  Peak files: $peak_count"
            log_message "  Size-selected BAMs: $bam_count"
            log_message "Skipping R analysis"
        else
            log_message "Found $peak_count peak files and $bam_count BAM files for analysis"

            # Validate sample metadata
            if ! validate_sample_metadata "$SPECIES_DIR"; then
                log_message "WARNING: Failed to validate sample metadata"
                if [[ "$R_ANALYSIS_REQUIRED" == "true" ]]; then
                    log_message "ERROR: Cannot proceed with R analysis"
                    exit 1
                else
                    log_message "Skipping R analysis"
                fi
            else
                # Run R analysis
                if run_r_analysis "$SPECIES_DIR"; then
                    create_checkpoint "r_analysis_complete"

                    # Send success email
                    if [[ "$SEND_EMAIL_UPDATES" == "true" ]] && [[ "$DRY_RUN" == "false" ]]; then
                        local motif_status=""
                        if [[ "$R_SKIP_MOTIFS" != "true" ]]; then
                            motif_status="- Motif analysis (HOMER)\n"
                        fi

                        send_email "R Analysis Complete: $SPECIES_NAME" \
                            "R analysis completed successfully for $SPECIES_NAME at $(date '+%Y-%m-%d %H:%M:%S')\n\nOutputs available in: $SPECIES_DIR/figures\n\nGenerated:\n- Venn diagrams\n- Correlation plots\n- Peak annotations\n- Genomic distribution\n- Feature tables\n- Consensus peaks\n${motif_status}\nPipeline Version: $PIPELINE_VERSION"
                    fi
                else
                    # Send failure email
                    if [[ "$SEND_EMAIL_UPDATES" == "true" ]] && [[ "$DRY_RUN" == "false" ]]; then
                        send_email "R Analysis Failed: $SPECIES_NAME" \
                            "R analysis failed for $SPECIES_NAME at $(date '+%Y-%m-%d %H:%M:%S')\n\nCheck error log: $PROCESSED_DIR/r_analysis_error.log\n\nPipeline continued successfully."
                    fi

                    if [[ "$R_ANALYSIS_REQUIRED" == "true" ]]; then
                        exit 1
                    fi
                fi
            fi
        fi
    fi
fi

##### Pipeline Completion #####

log_stage "Pipeline Summary for $SPECIES_NAME"

# Generate summary statistics
log_message "Generating summary statistics for $SPECIES_NAME..."
summary_file="$PROCESSED_DIR/pipeline_summary.txt"

if [[ "$DRY_RUN" == "false" ]]; then
    {
        echo "=== IISAGE ATAC-seq Pipeline Summary ==="
        echo "Pipeline Version: $PIPELINE_VERSION"
        echo "Species: $SPECIES_NAME"
        echo "Pipeline completed: $(date)"
        echo "==========================================="
        echo ""
        echo "Configuration:"
        echo "-------------"
        echo "Genome size: $GENOME_SIZE"
        echo "Reference: $REF_GENOME_FASTA"
        echo "Max insert size: $MAX_INSERT_SIZE"
        echo "Size selection cutoff: $SIZE_SELECTION_CUTOFF bp"
        echo "Minimum MAPQ: $MIN_MAPQ"
        echo "TSS analysis: ±${TSS_UPSTREAM}/${TSS_DOWNSTREAM}bp"
        echo "BigWig normalization: $BIGWIG_NORMALIZE"
        echo "BigWig bin size: ${BIGWIG_BINSIZE}bp"
        echo ""
        echo "Sample Summary:"
        echo "==============="
    } > "$summary_file"

    total_reads=0
    sample_count=0
    for bam in "$PROCESSED_DIR"/*_sorted.bam; do
        if [[ -f "$bam" ]]; then
            sample_name=$(basename "$bam" _sorted.bam)
            reads=$(samtools view -c "$bam")
            total_reads=$((total_reads + reads))
            sample_count=$((sample_count + 1))
            echo "$sample_name: $reads aligned reads" >> "$summary_file"
        fi
    done

    if [[ $sample_count -gt 0 ]]; then
        avg_reads=$((total_reads / sample_count))
        {
            echo ""
            echo "Total samples: $sample_count"
            echo "Total reads: $total_reads"
            echo "Average reads per sample: $avg_reads"
        } >> "$summary_file"
    fi

    if [[ "$RUN_PEAK_CALLING" == "true" ]]; then
        {
            echo ""
            echo "Peak Summary:"
            echo "============="
        } >> "$summary_file"

        total_peaks=0
        peak_samples=0
        for peak_file in "$PEAKS_DIR"/*_peaks.narrowPeak; do
            if [[ -f "$peak_file" ]]; then
                sample_name=$(basename "$peak_file" _peaks.narrowPeak)
                peak_count=$(wc -l < "$peak_file")
                total_peaks=$((total_peaks + peak_count))
                peak_samples=$((peak_samples + 1))
                echo "$sample_name: $peak_count peaks" >> "$summary_file"
            fi
        done

        if [[ $peak_samples -gt 0 ]]; then
            avg_peaks=$((total_peaks / peak_samples))
            {
                echo ""
                echo "Total peaks: $total_peaks"
                echo "Average peaks per sample: $avg_peaks"
            } >> "$summary_file"
        fi
    fi

    # R Analysis Summary
    if [[ "$RUN_R_ANALYSIS" == "true" ]] && check_checkpoint "r_analysis_complete"; then
        {
            echo ""
            echo "R Analysis Summary:"
            echo "==================="
        } >> "$summary_file"

        figures_dir="$SPECIES_DIR/figures"
        if [[ -d "$figures_dir" ]]; then
            figure_count=$(find "$figures_dir" -type f \( -name "*.png" -o -name "*.pdf" \) 2>/dev/null | wc -l)
            csv_count=$(find "$figures_dir" -type f -name "*.csv" 2>/dev/null | wc -l)

            echo "Figures generated: $figure_count" >> "$summary_file"
            echo "Summary tables: $csv_count" >> "$summary_file"
            echo "Output directory: $figures_dir" >> "$summary_file"

            if [[ "$R_SKIP_MOTIFS" != "true" ]]; then
                motif_dirs=$(find "$figures_dir" -type d -name "*_motifs" 2>/dev/null | wc -l)
                if [[ $motif_dirs -gt 0 ]]; then
                    echo "Motif analysis: $motif_dirs samples" >> "$summary_file"
                fi
            fi
        fi
    fi

    {
        echo ""
        echo "Software versions logged in: $VERSION_FILE"
    } >> "$summary_file"
fi

log_message "ATAC-seq pipeline completed successfully for $SPECIES_NAME!"
log_message "Output directories:"
log_message "  - Processed BAMs: $PROCESSED_DIR"
log_message "  - Peak calling results: $PEAKS_DIR"
log_message "  - DeepTools visualizations: $DEEPTOOLS_DIR"
if [[ "$RUN_R_ANALYSIS" == "true" ]] && check_checkpoint "r_analysis_complete"; then
    log_message "  - R analysis figures: $SPECIES_DIR/figures"
fi
log_message "Summary statistics saved to: $summary_file"

# Send completion email
if [[ "$SEND_EMAIL_UPDATES" == "true" ]] && [[ "$DRY_RUN" == "false" ]]; then
    summary_content=$(cat "$summary_file")

    # Build output directories message
    output_dirs="Output directories:\n- Processed BAMs: $PROCESSED_DIR\n- Peak calling results: $PEAKS_DIR\n- DeepTools visualizations: $DEEPTOOLS_DIR"
    if [[ "$RUN_R_ANALYSIS" == "true" ]] && check_checkpoint "r_analysis_complete"; then
        output_dirs="${output_dirs}\n- R analysis figures: $SPECIES_DIR/figures"
    fi

    send_email "ATAC-seq Pipeline Completed: $SPECIES_NAME" \
               "Pipeline v$PIPELINE_VERSION completed for $SPECIES_NAME at $(date '+%Y-%m-%d %H:%M:%S')\n\n$summary_content\n\n${output_dirs}"
fi

log_message "Pipeline execution completed at $(date)"

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
#SBATCH --mail-user=decenalo@msu.edu
#SBATCH --error=%x-%j.SLURMerr

####################################################################################
# Script: multi_species_atac_pipeline.sh
# Purpose: Multi-species ATAC-seq pipeline with YAML configuration
# Author: Louis Paul Decena-Segarra, PhD
# Created: 2025-06-02
# Usage:
#   Single species: sbatch --export=SPECIES_NAME=species_name multi_species_atac_pipeline.sh
#   All species: sbatch --export=DEPLOY_ALL=true multi_species_atac_pipeline.sh
# Pipeline stages: 1) Quality control & trimming, 2) Alignment & BAM processing,
#                  3) Peak calling with size selection, 4) DeepTools visualization
#
# Directory structure expected:
# $BASE_DIR/
# ├── species1/
# │   ├── config.yaml
# │   ├── samples/
# │   │   ├── sample1/
# │   │   └── sample2/
# │   └── reference/
# └── species2/
#     ├── config.yaml
#     ├── samples/
#     └── reference/
####################################################################################

set -euo pipefail

##### Global Configuration #####

# Base directory containing all species folders. Full path prefered
BASE_DIR="/mnt/scratch/decenalo/atac_seq_projects"

# Processing parameters (can be overridden by YAML)
DEFAULT_THREADS=16
DEFAULT_MAX_INSERT_SIZE=2000
DEFAULT_SIZE_SELECTION_CUTOFF=100
DEFAULT_MIN_MAPQ=5
DEFAULT_TRIM_SETTINGS="ILLUMINACLIP:{ADAPTERS_PATH}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:25"
DEFAULT_BOWTIE2_OPTS="--very-sensitive --no-mixed --no-discordant --dovetail"


# Pipeline control
# Changing these parameters from true/false will change which steps are performed
RUN_TRIMMING=true
RUN_ALIGNMENT=true
RUN_PEAK_CALLING=true
RUN_DEEPTOOLS=true
CLEANUP_INTERMEDIATE=false

# Email notifications
#Change this to recieve updates to your email (please, don't keep my email in there)
SEND_EMAIL_UPDATES=true
EMAIL_ADDRESS="decenalo@msu.edu"

##### Functions #####
#These functions do different things for the pipeline

#Create log messages
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log_stage() {
    echo ""
    echo "==========================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] STAGE: $1"
    echo "==========================================="

    if [[ "$SEND_EMAIL_UPDATES" == "true" ]]; then
        send_email "ATAC-seq Pipeline ($CURRENT_SPECIES): $1 Started" \
                   "Stage '$1' has begun for species $CURRENT_SPECIES at $(date '+%Y-%m-%d %H:%M:%S')\nJob ID: $SLURM_JOB_ID\nNode: $SLURM_NODELIST"
    fi
}

#Send emails
send_email() {
    if [[ "$SEND_EMAIL_UPDATES" == "true" ]] && command -v mail >/dev/null 2>&1; then
        local subject="$1"
        local message="$2"
        echo -e "$message" | mail -s "$subject" "$EMAIL_ADDRESS" 2>/dev/null || true
    fi
}

#Check if files exists
check_file_exists() {
    if [[ ! -f "$1" ]]; then
        log_message "ERROR: Required file not found: $1"
        exit 1
    fi
}

#Check if directories exist
check_dir_exists() {
    if [[ ! -d "$1" ]]; then
        log_message "ERROR: Required directory not found: $1"
        exit 1
    fi
}

#Parse the yaml (configuration file) for each species
parse_yaml() {
    local yaml_file="$1"
    local prefix="$2"

    # Simple YAML parser using awk
    awk -F': ' -v prefix="$prefix" '
    /^[[:space:]]*#/ { next }
    /^[[:space:]]*$/ { next }
    /^[^[:space:]]/ {
        gsub(/[[:space:]]*$/, "", $1)
        if (NF >= 2) {
            gsub(/[[:space:]]*/, "", $2)
            gsub(/["'"'"']/, "", $2)
            print prefix $1 "=" $2
        }
    }
    ' "$yaml_file"
}

#Load the species configuration from the yaml file
load_species_config() {
    local species_dir="$1"
    local config_file="$species_dir/config.yaml"

    if [[ ! -f "$config_file" ]]; then
        log_message "ERROR: Config file not found: $config_file"
        return 1
    fi

    log_message "Loading configuration from: $config_file"

    # Parse YAML and set variables
    while IFS= read -r line; do
        if [[ "$line" =~ ^[A-Z_]+=.* ]]; then
            eval "export $line"
            log_message "Config: $line"
        fi
    done < <(parse_yaml "$config_file" "")

    # Set species-specific paths
    export SPECIES_DIR="$species_dir"
    export SAMPLES_DIR="$species_dir/samples"
    export PROCESSED_DIR="$species_dir/processed_files"
    export PEAKS_DIR="$species_dir/peaks_output"
    export DEEPTOOLS_DIR="$species_dir/deeptools_output"

    # Validate required config variables
    local required_vars=("GENOME_SIZE" "REF_GENOME_FASTA" "SPECIES_NAME")
    for var in "${required_vars[@]}"; do
        if [[ -z "${!var:-}" ]]; then
            log_message "ERROR: Required variable $var not found in config.yaml"
            return 1
        fi
    done

    # Set optional variables with defaults
    export THREADS="${THREADS:-$DEFAULT_THREADS}"
    export MAX_INSERT_SIZE="${MAX_INSERT_SIZE:-$DEFAULT_MAX_INSERT_SIZE}"
    export SIZE_SELECTION_CUTOFF="${SIZE_SELECTION_CUTOFF:-$DEFAULT_SIZE_SELECTION_CUTOFF}"
    export MIN_MAPQ="${MIN_MAPQ:-$DEFAULT_MIN_MAPQ}"
    export BOWTIE2_OPTS="${BOWTIE2_OPTS:-$DEFAULT_BOWTIE2_OPTS}"

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

#Deploy species job (submiting and individual SLURM job for each species)
deploy_species_job() {
    local species_name="$1"
    local species_dir="$BASE_DIR/$species_name"

    if [[ ! -d "$species_dir" ]]; then
        log_message "WARNING: Species directory not found: $species_dir"
        return 1
    fi

    log_message "Deploying job for species: $species_name"

    # Submit job for this species
    sbatch --export=SPECIES_NAME="$species_name" \
           --job-name="atac-${species_name}" \
           --output="atac-${species_name}-%j.SLURMout" \
           --error="atac-${species_name}-%j.SLURMerr" \
           "$0"

    log_message "Job submitted for $species_name"
}

#Remove inermediate files
cleanup_files() {
    local files=("$@")
    if [[ "$CLEANUP_INTERMEDIATE" == "true" ]]; then
        for file in "${files[@]}"; do
            if [[ -f "$file" ]]; then
                log_message "Cleaning up: $(basename "$file")"
                rm "$file"
            fi
        done
    fi
}

##### Main Pipeline Logic #####

# Change to working directory
cd /mnt/scratch/decenalo/

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
               "Deployed $species_count species jobs at $(date '+%Y-%m-%d %H:%M:%S')"
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

# Load modules
log_message "Loading required modules"
module purge
module load Trimmomatic Bowtie2 SAMtools MACS2/2.2.9.1-foss-2022b deepTools

# Validate paths
check_dir_exists "$SAMPLES_DIR"
check_file_exists "$REF_GENOME_FASTA"

# Create output directories
mkdir -p "$PROCESSED_DIR" "$PEAKS_DIR" "$DEEPTOOLS_DIR"

# Log configuration summary
log_message "Configuration summary for $SPECIES_NAME:"
log_message "  Genome size: $GENOME_SIZE"
log_message "  Reference: $REF_GENOME_FASTA"
log_message "  Samples directory: $SAMPLES_DIR"
log_message "  Threads: $THREADS"

##### STAGE 1: Trimming & Alignment #####

if [[ "$RUN_TRIMMING" == "true" ]] || [[ "$RUN_ALIGNMENT" == "true" ]]; then
    log_stage "Quality Control, Trimming & Alignment"

    # Check/build Bowtie2 index
    if [[ ! -f "${REF_GENOME}.1.bt2" ]]; then
        log_message "Building Bowtie2 index for $SPECIES_NAME..."
        mkdir -p "$(dirname "$REF_GENOME")"
        bowtie2-build "$REF_GENOME_FASTA" "$REF_GENOME"
        log_message "Index created successfully"
    else
        log_message "Bowtie2 index found for $SPECIES_NAME"
    fi

    # Count and process samples
    TOTAL_SAMPLES=$(find "$SAMPLES_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
    log_message "Found $TOTAL_SAMPLES samples for $SPECIES_NAME"

    CURRENT_SAMPLE=0
    for SAMPLE_DIR in "$SAMPLES_DIR"/*; do
        if [[ ! -d "$SAMPLE_DIR" ]]; then
            continue
        fi

        SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        ((CURRENT_SAMPLE++))

        log_message "Processing sample $CURRENT_SAMPLE/$TOTAL_SAMPLES: $SAMPLE_NAME ($SPECIES_NAME)"

        # Check for input files
        R1_FILES=("$SAMPLE_DIR"/*_1.fq.gz)
        R2_FILES=("$SAMPLE_DIR"/*_2.fq.gz)

        if [[ ! -e "${R1_FILES[0]}" ]] || [[ ! -e "${R2_FILES[0]}" ]]; then
            log_message "WARNING: No paired FASTQ files found in $SAMPLE_DIR, skipping"
            continue
        fi

        # Define file paths
        BAM_FILE="$PROCESSED_DIR/${SAMPLE_NAME}_sorted.bam"

        # Skip if BAM already exists
        if [[ -f "$BAM_FILE" ]] && [[ -f "${BAM_FILE}.bai" ]] && [[ "$RUN_ALIGNMENT" == "false" ]]; then
            log_message "BAM file exists for $SAMPLE_NAME, skipping alignment"
            continue
        fi

        if [[ "$RUN_TRIMMING" == "true" ]] && [[ "$RUN_ALIGNMENT" == "true" ]]; then
            # File paths
            MERGED_R1="$PROCESSED_DIR/${SAMPLE_NAME}_R1_merged.fastq"
            MERGED_R2="$PROCESSED_DIR/${SAMPLE_NAME}_R2_merged.fastq"
            PAIRED_R1="$PROCESSED_DIR/${SAMPLE_NAME}_R1_paired.fastq"
            PAIRED_R2="$PROCESSED_DIR/${SAMPLE_NAME}_R2_paired.fastq"
            UNPAIRED_R1="$PROCESSED_DIR/${SAMPLE_NAME}_R1_unpaired.fastq"
            UNPAIRED_R2="$PROCESSED_DIR/${SAMPLE_NAME}_R2_unpaired.fastq"
            SAM_FILE="$PROCESSED_DIR/${SAMPLE_NAME}_aligned.sam"

            # Merge FASTQ files
            log_message "Merging FASTQ files for $SAMPLE_NAME"
            zcat "$SAMPLE_DIR"/*_1.fq.gz > "$MERGED_R1"
            zcat "$SAMPLE_DIR"/*_2.fq.gz > "$MERGED_R2"

            # Trim adapters
            log_message "Trimming adapters for $SAMPLE_NAME"
            java -jar "$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar" PE -threads "$THREADS" \
                "$MERGED_R1" "$MERGED_R2" \
                "$PAIRED_R1" "$UNPAIRED_R1" "$PAIRED_R2" "$UNPAIRED_R2" \
                $TRIM_SETTINGS

            cleanup_files "$MERGED_R1" "$MERGED_R2"

            # Alignment
            log_message "Aligning reads for $SAMPLE_NAME"
            bowtie2 -x "$REF_GENOME" -1 "$PAIRED_R1" -2 "$PAIRED_R2" \
                    -S "$SAM_FILE" $BOWTIE2_OPTS \
                    -X "$MAX_INSERT_SIZE" -k 2 -p "$THREADS"

            # Convert to BAM, sort and index
            log_message "Converting to BAM and sorting for $SAMPLE_NAME"
            samtools view -bS "$SAM_FILE" | samtools sort -@ "$THREADS" -o "$BAM_FILE"
            samtools index "$BAM_FILE"

            cleanup_files "$PAIRED_R1" "$PAIRED_R2" "$UNPAIRED_R1" "$UNPAIRED_R2" "$SAM_FILE"
        fi
    done
fi

##### STAGE 2: Peak Calling with Size Selection #####

if [[ "$RUN_PEAK_CALLING" == "true" ]]; then
    log_stage "Peak Calling with Size Selection"

    for BAM_FILE in "$PROCESSED_DIR"/*_sorted.bam; do
        if [[ ! -f "$BAM_FILE" ]]; then
            log_message "No BAM files found in $PROCESSED_DIR"
            break
        fi

        SAMPLE_NAME=$(basename "$BAM_FILE" _sorted.bam)
        log_message "Processing peaks for $SAMPLE_NAME ($SPECIES_NAME)"

        # Check if peaks already exist
        if [[ -f "$PEAKS_DIR/${SAMPLE_NAME}_peaks.narrowPeak" ]]; then
            log_message "Peaks already exist for $SAMPLE_NAME, skipping"
            continue
        fi

        # Size selection and filtering
        log_message "Performing size selection and filtering for $SAMPLE_NAME"
        SIZE_SELECTED_BAM="${PEAKS_DIR}/${SAMPLE_NAME}_size_selected_${SIZE_SELECTION_CUTOFF}.bam"

        samtools view -@ "$THREADS" -h "$BAM_FILE" | \
        awk -v LEN="$SIZE_SELECTION_CUTOFF" '{
            if ($9 <= LEN && $9 >= -LEN && $9 != 0 || $1 ~ /^@/) print $0
        }' | \
        grep -v "_random\|chrUn\|chrM" | \
        samtools view -@ "$THREADS" -uh -q "$MIN_MAPQ" - | \
        samtools sort -@ "$THREADS" -o "$SIZE_SELECTED_BAM"

        samtools index -@ "$THREADS" "$SIZE_SELECTED_BAM"

        # Call peaks with MACS2
        log_message "Calling peaks for $SAMPLE_NAME"
        macs2 callpeak -t "$SIZE_SELECTED_BAM" \
                       -f BAM \
                       -g "$GENOME_SIZE" \
                       --nomodel \
                       --nolambda \
                       --keep-dup all \
                       --call-summits \
                       -n "$SAMPLE_NAME" \
                       --outdir "$PEAKS_DIR"

        log_message "Peaks called for $SAMPLE_NAME"
    done
fi

##### STAGE 3: DeepTools Visualization #####

if [[ "$RUN_DEEPTOOLS" == "true" ]]; then
    log_stage "DeepTools Visualization"

    # Collect all size-selected BAM files
    SIZE_SELECTED_BAMS=("$PEAKS_DIR"/*_size_selected_${SIZE_SELECTION_CUTOFF}.bam)

    if [[ ${#SIZE_SELECTED_BAMS[@]} -eq 0 ]] || [[ ! -f "${SIZE_SELECTED_BAMS[0]}" ]]; then
        log_message "No size-selected BAM files found for DeepTools analysis"
    else
        log_message "Running DeepTools analysis on ${#SIZE_SELECTED_BAMS[@]} samples for $SPECIES_NAME"

        # Create sample labels
        SAMPLE_LABELS=""
        for bam in "${SIZE_SELECTED_BAMS[@]}"; do
            sample_name=$(basename "$bam" "_size_selected_${SIZE_SELECTION_CUTOFF}.bam")
            SAMPLE_LABELS="$SAMPLE_LABELS $sample_name"
        done

        # Look for TSS file in species directory
        TSS_BED="$SPECIES_DIR/reference/tss_regions.bed"

        if [[ -f "$TSS_BED" ]]; then
            log_message "Computing matrix around TSS regions"
            computeMatrix reference-point \
                --referencePoint TSS \
                -b 2000 -a 2000 \
                -R "$TSS_BED" \
                -S "${SIZE_SELECTED_BAMS[@]}" \
                --skipZeros \
                -o "$DEEPTOOLS_DIR/matrix_TSS.gz" \
                --samplesLabel $SAMPLE_LABELS \
                -p "$THREADS"

            # Create heatmap
            log_message "Creating TSS heatmap"
            plotHeatmap -m "$DEEPTOOLS_DIR/matrix_TSS.gz" \
                       -out "$DEEPTOOLS_DIR/TSS_heatmap.png" \
                       --colorMap viridis \
                       --whatToShow 'heatmap and colorbar'
        else
            log_message "TSS BED file not found at $TSS_BED, skipping TSS analysis"
        fi

        # Compute coverage correlation
        log_message "Computing sample correlation"
        multiBamSummary bins \
            --bamfiles "${SIZE_SELECTED_BAMS[@]}" \
            --labels $SAMPLE_LABELS \
            -out "$DEEPTOOLS_DIR/multiBamSummary.npz" \
            -p "$THREADS"

        # Plot correlation
        plotCorrelation \
            -in "$DEEPTOOLS_DIR/multiBamSummary.npz" \
            --corMethod pearson \
            --skipZeros \
            --plotTitle "Sample Correlation - $SPECIES_NAME" \
            --whatToPlot heatmap \
            --colorMap RdYlBu \
            --plotNumbers \
            -o "$DEEPTOOLS_DIR/correlation_heatmap.png"

        # Create bigWig files
        log_message "Creating bigWig files for genome browser visualization"
        for bam in "${SIZE_SELECTED_BAMS[@]}"; do
            sample_name=$(basename "$bam" "_size_selected_${SIZE_SELECTION_CUTOFF}.bam")
            log_message "Creating bigWig for $sample_name"
            bamCoverage -b "$bam" \
                       -o "$DEEPTOOLS_DIR/${sample_name}.bw" \
                       --normalizeUsing RPKM \
                       --binSize 10 \
                       -p "$THREADS"
        done
    fi
fi

##### Pipeline Completion #####

log_stage "Pipeline Summary for $SPECIES_NAME"

# Generate summary statistics
log_message "Generating summary statistics for $SPECIES_NAME..."
SUMMARY_FILE="$PROCESSED_DIR/pipeline_summary.txt"
echo "Species: $SPECIES_NAME" > "$SUMMARY_FILE"
echo "Pipeline completed: $(date)" >> "$SUMMARY_FILE"
echo "===========================================" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
echo "Sample Summary:" >> "$SUMMARY_FILE"
echo "===============" >> "$SUMMARY_FILE"

total_reads=0
sample_count=0
for bam in "$PROCESSED_DIR"/*_sorted.bam; do
    if [[ -f "$bam" ]]; then
        sample_name=$(basename "$bam" _sorted.bam)
        reads=$(samtools view -c "$bam")
        total_reads=$((total_reads + reads))
        sample_count=$((sample_count + 1))
        echo "$sample_name: $reads aligned reads" >> "$SUMMARY_FILE"
    fi
done

if [[ $sample_count -gt 0 ]]; then
    avg_reads=$((total_reads / sample_count))
    echo "" >> "$SUMMARY_FILE"
    echo "Total samples: $sample_count" >> "$SUMMARY_FILE"
    echo "Total reads: $total_reads" >> "$SUMMARY_FILE"
    echo "Average reads per sample: $avg_reads" >> "$SUMMARY_FILE"
fi

if [[ "$RUN_PEAK_CALLING" == "true" ]]; then
    echo "" >> "$SUMMARY_FILE"
    echo "Peak Summary:" >> "$SUMMARY_FILE"
    echo "=============" >> "$SUMMARY_FILE"
    total_peaks=0
    peak_samples=0
    for peak_file in "$PEAKS_DIR"/*_peaks.narrowPeak; do
        if [[ -f "$peak_file" ]]; then
            sample_name=$(basename "$peak_file" _peaks.narrowPeak)
            peak_count=$(wc -l < "$peak_file")
            total_peaks=$((total_peaks + peak_count))
            peak_samples=$((peak_samples + 1))
            echo "$sample_name: $peak_count peaks" >> "$SUMMARY_FILE"
        fi
    done

    if [[ $peak_samples -gt 0 ]]; then
        avg_peaks=$((total_peaks / peak_samples))
        echo "" >> "$SUMMARY_FILE"
        echo "Total peaks: $total_peaks" >> "$SUMMARY_FILE"
        echo "Average peaks per sample: $avg_peaks" >> "$SUMMARY_FILE"
    fi
fi

log_message "ATAC-seq pipeline completed successfully for $SPECIES_NAME!"
log_message "Output directories:"
log_message "  - Processed BAMs: $PROCESSED_DIR"
log_message "  - Peak calling results: $PEAKS_DIR"
log_message "  - DeepTools visualizations: $DEEPTOOLS_DIR"
log_message "Summary statistics saved to: $SUMMARY_FILE"

# Send completion email
if [[ "$SEND_EMAIL_UPDATES" == "true" ]]; then
    SUMMARY_CONTENT=$(cat "$SUMMARY_FILE")
    send_email "ATAC-seq Pipeline Completed: $SPECIES_NAME" \
               "Pipeline completed for $SPECIES_NAME at $(date '+%Y-%m-%d %H:%M:%S')\n\n$SUMMARY_CONTENT\n\nOutput directories:\n- Processed BAMs: $PROCESSED_DIR\n- Peak calling results: $PEAKS_DIR\n- DeepTools visualizations: $DEEPTOOLS_DIR"
fi

log_message "Pipeline execution completed at $(date)"

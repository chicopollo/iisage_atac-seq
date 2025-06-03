#!/bin/bash --login

##### SLURM Resource Requests #####

#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G

##### SLURM Administrative Settings #####

#SBATCH --job-name=2025_Cpicta_atac-seq_pipeline
#SBATCH --output=%x-%j.SLURMout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=x_user@university.edu
#SBATCH --error=%x-%j.SLURMerr

####################################################################################
# Script: unified_atac_seq_pipeline.sh
# Purpose: Complete ATAC-seq pipeline - trimming, alignment, peak calling, and visualization
# Author: x_user, PhD
# Created: 2025-06-02
# Usage: sbatch unified_atac_seq_pipeline.sh
# Pipeline stages: 1) Quality control & trimming, 2) Alignment & BAM processing,
#                  3) Peak calling with size selection, 4) DeepTools visualization
# Notes: Customize variables section before running!!
###################################################################################

set -euo pipefail  # Exit on error, undefined vars, pipe failures

##### Configuration Variables #####

# Main directories - UPDATE THESE FOR YOUR SETUP
SAMPLES_DIR="/mnt/scratch/x_user/atac_seq"
PROCESSED_DIR="/mnt/scratch/x_user/processed_files"
PEAKS_DIR="/mnt/scratch/x_user/peaks_output"
DEEPTOOLS_DIR="/mnt/scratch/x_user/deeptools_output"
REF_GENOME="/mnt/scratch/x_user/genome_index_c_picta"
REF_GENOME_FASTA="/mnt/scratch/x_user/cpicta_assembly/GCF_011386835.1_ASM1138683v2_genomic.fna"
ADAPTERS_PATH="/mnt/home/x_user/Documents/adapters/novogene.fa"

# Species-specific parameters
GENOME_SIZE="2.4e9"  # Chrysemys picta genome size for MACS2

# Processing parameters
THREADS=16
MAX_INSERT_SIZE=2000
SIZE_SELECTION_CUTOFF=100  # Fragment size cutoff for ATAC-seq
MIN_MAPQ=5  # Minimum mapping quality

# Trimming settings
TRIM_SETTINGS="ILLUMINACLIP:${ADAPTERS_PATH}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:25"

# Bowtie2 alignment options
BOWTIE2_OPTS="--very-sensitive --no-mixed --no-discordant --dovetail"

# Pipeline control - set to false to skip stages
RUN_TRIMMING=true
RUN_ALIGNMENT=true
RUN_PEAK_CALLING=true
RUN_DEEPTOOLS=true

# Cleanup options
CLEANUP_INTERMEDIATE=false

# Email notifications
SEND_EMAIL_UPDATES=true
EMAIL_ADDRESS="x_user@university.edu"

##### Functions #####

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log_stage() {
    echo ""
    echo "==========================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] STAGE: $1"
    echo "==========================================="

    # Send email notification for major stages
    if [[ "$SEND_EMAIL_UPDATES" == "true" ]]; then
        send_email "ATAC-seq Pipeline: $1 Started" "Stage '$1' has begun at $(date '+%Y-%m-%d %H:%M:%S')\nJob ID: $SLURM_JOB_ID\nNode: $SLURM_NODEID"
    fi
}

send_email() {
    if [[ "$SEND_EMAIL_UPDATES" == "true" ]] && command -v mail >/dev/null 2>&1; then
        local subject="$1"
        local message="$2"
        echo -e "$message" | mail -s "$subject" "$EMAIL_ADDRESS" 2>/dev/null || true
    fi
}

check_file_exists() {
    if [[ ! -f "$1" ]]; then
        log_message "ERROR: Required file not found: $1"
        exit 1
    fi
}

check_dir_exists() {
    if [[ ! -d "$1" ]]; then
        log_message "ERROR: Required directory not found: $1"
        exit 1
    fi
}

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

##### Main Pipeline #####

log_stage "ATAC-seq Pipeline Initialization"
log_message "Pipeline stages enabled:"
log_message "  - Trimming & QC: $RUN_TRIMMING"
log_message "  - Alignment: $RUN_ALIGNMENT"
log_message "  - Peak calling: $RUN_PEAK_CALLING"
log_message "  - DeepTools visualization: $RUN_DEEPTOOLS"

# Load all required modules
log_message "Loading required modules"
module purge
module load Trimmomatic Bowtie2 SAMtools MACS2/2.2.9.1-foss-2022b deepTools

# Change to working directory
cd /mnt/scratch/x_user/

# Validate inputs
log_message "Validating input paths"
check_dir_exists "$SAMPLES_DIR"
check_file_exists "$REF_GENOME_FASTA"
check_file_exists "$ADAPTERS_PATH"

# Create output directories
mkdir -p "$PROCESSED_DIR" "$PEAKS_DIR" "$DEEPTOOLS_DIR"

# Count samples
TOTAL_SAMPLES=$(find "$SAMPLES_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
log_message "Found $TOTAL_SAMPLES samples to process"

##### STAGE 1: Trimming & Alignment #####

if [[ "$RUN_TRIMMING" == "true" ]] || [[ "$RUN_ALIGNMENT" == "true" ]]; then
    log_stage "Quality Control, Trimming & Alignment"

    # Check/build Bowtie2 index
    if [[ ! -f "${REF_GENOME}.1.bt2" ]]; then
        log_message "Building Bowtie2 index..."
        bowtie2-build "$REF_GENOME_FASTA" "$REF_GENOME"
        log_message "Index created successfully"
    else
        log_message "Bowtie2 index found"
    fi

    CURRENT_SAMPLE=0
    for SAMPLE_DIR in "$SAMPLES_DIR"/*; do
        if [[ ! -d "$SAMPLE_DIR" ]]; then
            continue
        fi

        SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        ((CURRENT_SAMPLE++))

        log_message "Processing sample $CURRENT_SAMPLE/$TOTAL_SAMPLES: $SAMPLE_NAME"

        # Check for input files
        R1_FILES=("$SAMPLE_DIR"/*_1.fq.gz)
        R2_FILES=("$SAMPLE_DIR"/*_2.fq.gz)

        if [[ ! -e "${R1_FILES[0]}" ]] || [[ ! -e "${R2_FILES[0]}" ]]; then
            log_message "WARNING: No paired FASTQ files found in $SAMPLE_DIR, skipping"
            continue
        fi

        # Define file paths
        BAM_FILE="$PROCESSED_DIR/${SAMPLE_NAME}_sorted.bam"

        # Skip if BAM already exists and we're not rerunning alignment
        if [[ -f "$BAM_FILE" ]] && [[ -f "${BAM_FILE}.bai" ]] && [[ "$RUN_ALIGNMENT" == "false" ]]; then
            log_message "BAM file exists for $SAMPLE_NAME, skipping alignment"
            continue
        fi

        if [[ "$RUN_TRIMMING" == "true" ]] && [[ "$RUN_ALIGNMENT" == "true" ]]; then
            # File paths for trimming
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
        log_message "Processing peaks for $SAMPLE_NAME"

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
        log_message "Running DeepTools analysis on ${#SIZE_SELECTED_BAMS[@]} samples"

        # Create sample labels
        SAMPLE_LABELS=""
        for bam in "${SIZE_SELECTED_BAMS[@]}"; do
            sample_name=$(basename "$bam" "_size_selected_${SIZE_SELECTION_CUTOFF}.bam")
            SAMPLE_LABELS="$SAMPLE_LABELS $sample_name"
        done

        # Compute matrix for TSS regions
        log_message "Computing matrix around TSS regions"
        # Note: You'll need a BED file with TSS coordinates for your species
        # This is a placeholder - adjust according to your annotation files
        TSS_BED="$PROCESSED_DIR/tss_regions.bed"  # You need to provide this

        if [[ -f "$TSS_BED" ]]; then
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
            log_message "WARNING: TSS BED file not found at $TSS_BED, skipping TSS analysis"
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
            --plotTitle "Sample Correlation" \
            --whatToPlot heatmap \
            --colorMap RdYlBu \
            --plotNumbers \
            -o "$DEEPTOOLS_DIR/correlation_heatmap.png"

        # Create bigWig files for visualization
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

log_stage "Pipeline Summary"
log_message "ATAC-seq pipeline completed successfully!"
log_message "Output directories:"
log_message "  - Processed BAMs: $PROCESSED_DIR"
log_message "  - Peak calling results: $PEAKS_DIR"
log_message "  - DeepTools visualizations: $DEEPTOOLS_DIR"

# Generate summary statistics
log_message "Generating summary statistics..."
echo "Sample Summary:" > "$PROCESSED_DIR/pipeline_summary.txt"
echo "===============" >> "$PROCESSED_DIR/pipeline_summary.txt"
for bam in "$PROCESSED_DIR"/*_sorted.bam; do
    if [[ -f "$bam" ]]; then
        sample_name=$(basename "$bam" _sorted.bam)
        total_reads=$(samtools view -c "$bam")
        echo "$sample_name: $total_reads aligned reads" >> "$PROCESSED_DIR/pipeline_summary.txt"
    fi
done

if [[ "$RUN_PEAK_CALLING" == "true" ]]; then
    echo "" >> "$PROCESSED_DIR/pipeline_summary.txt"
    echo "Peak Summary:" >> "$PROCESSED_DIR/pipeline_summary.txt"
    echo "=============" >> "$PROCESSED_DIR/pipeline_summary.txt"
    for peak_file in "$PEAKS_DIR"/*_peaks.narrowPeak; do
        if [[ -f "$peak_file" ]]; then
            sample_name=$(basename "$peak_file" _peaks.narrowPeak)
            peak_count=$(wc -l < "$peak_file")
            echo "$sample_name: $peak_count peaks" >> "$PROCESSED_DIR/pipeline_summary.txt"
        fi
    done
fi

log_message "Summary statistics saved to: $PROCESSED_DIR/pipeline_summary.txt"
log_message "Pipeline execution completed at $(date)"

# Send completion email with summary
if [[ "$SEND_EMAIL_UPDATES" == "true" ]]; then
    SUMMARY_CONTENT=$(cat "$PROCESSED_DIR/pipeline_summary.txt")
    send_email "ATAC-seq Pipeline Completed Successfully" "Pipeline completed at $(date '+%Y-%m-%d %H:%M:%S')\n\nSummary:\n$SUMMARY_CONTENT\n\nOutput directories:\n- Processed BAMs: $PROCESSED_DIR\n- Peak calling results: $PEAKS_DIR\n- DeepTools visualizations: $DEEPTOOLS_DIR"
fi

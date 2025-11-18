#!/bin/bash

#######################################################################################
# Script: atac_seq_pipeline_deployment.sh
# Purpose: Helper script to deploy ATAC-seq pipeline jobs with validation
# Author: Louis Paul Decena-Segarra, PhD
# Version: 2.0.0
# Created: 2025-06-02
# Modified: 2025-11-17
#
# Two-Job Workflow:
#   Job 1 (Array): Processes samples in parallel (trim, align, create BAMs)
#   Job 2 (Single): Runs downstream analysis (peak calling, deepTools, R analysis)
#                   Automatically submitted with dependency on Job 1 completion
#
# Usage Examples:
#   ./atac_seq_pipeline_deployment.sh --species Chrysemys_picta      # Deploy single species
#   ./atac_seq_pipeline_deployment.sh --deploy-all                   # Deploy all valid species
#   ./atac_seq_pipeline_deployment.sh --validate --species Chrysemys_picta # Validate first
#   ./atac_seq_pipeline_deployment.sh --validate                     # Validate all
#######################################################################################

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="/mnt/scratch/decenalo/atac_seq_projects"
PIPELINE_SCRIPT="$SCRIPT_DIR/multi_species_atac_pipeline_v2.5.sh"
R_MODULE="R/4.2.2-foss-2022b"
DEFAULT_CPUS=4
DEFAULT_MEM_PER_CPU=4  # GB

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --species NAME       Deploy pipeline for specific species"
    echo "  --deploy-all         Deploy pipeline for all valid species"
    echo "  --validate           Validate configuration before deployment"
    echo "  --cpus N             CPUs per task (default: $DEFAULT_CPUS)"
    echo "  --mem N              Memory per CPU in GB (default: ${DEFAULT_MEM_PER_CPU}G)"
    echo "  --help               Show this help message"
    echo ""
    echo "Examples:"
    echo "  # Validate configuration first"
    echo "  $0 --validate --species Chrysemys_picta"
    echo ""
    echo "  # Deploy single species with default resources"
    echo "  $0 --species Chrysemys_picta"
    echo ""
    echo "  # Deploy with custom resources"
    echo "  $0 --species Chrysemys_picta --cpus 8 --mem 8"
    echo ""
    echo "  # Validate all species"
    echo "  $0 --validate"
    echo ""
    echo "  # Deploy all valid species"
    echo "  $0 --deploy-all"
}

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[✓]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[✗]${NC} $1"
}

log_section() {
    echo ""
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}$1${NC}"
    echo -e "${CYAN}========================================${NC}"
}

check_prerequisites() {
    # Check if base directory exists
    if [[ ! -d "$BASE_DIR" ]]; then
        log_error "Base directory not found: $BASE_DIR"
        return 1
    fi

    # Check if pipeline script exists
    if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
        log_error "Pipeline script not found: $PIPELINE_SCRIPT"
        return 1
    fi

    # Check if SLURM is available
    if ! command -v sbatch >/dev/null 2>&1; then
        log_error "SLURM (sbatch) not available"
        return 1
    fi

    return 0
}

validate_species_config() {
    local species=$1
    local species_dir="$BASE_DIR/$species"
    local config_file="$species_dir/config.yaml"
    local errors=0

    log_info "Validating configuration for: $species"

    # Check species directory exists
    if [[ ! -d "$species_dir" ]]; then
        log_error "Species directory not found: $species_dir"
        return 1
    fi

    # Check config.yaml exists
    if [[ ! -f "$config_file" ]]; then
        log_error "config.yaml not found: $config_file"
        return 1
    fi

    # Check required fields in config.yaml
    local required_fields=("SPECIES_NAME" "GENOME_SIZE" "REF_GENOME_FASTA")
    for field in "${required_fields[@]}"; do
        if ! grep -q "^${field}:" "$config_file"; then
            log_error "Missing required field in config.yaml: $field"
            ((errors++))
        fi
    done

    # Check genome size format
    local genome_size=$(grep "^GENOME_SIZE:" "$config_file" | awk '{print $2}')
    if [[ -n "$genome_size" ]]; then
        if [[ ! "$genome_size" =~ ^[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]+)?$ ]]; then
            log_error "Invalid genome size format: $genome_size"
            ((errors++))
        fi
    fi

    # Check reference FASTA exists
    local ref_fasta=$(grep "^REF_GENOME_FASTA:" "$config_file" | cut -d: -f2- | tr -d ' "')
    if [[ -n "$ref_fasta" ]]; then
        # Handle relative paths
        if [[ ! "$ref_fasta" = /* ]]; then
            ref_fasta="$species_dir/$ref_fasta"
        fi
        if [[ ! -f "$ref_fasta" ]]; then
            log_warning "Reference FASTA not found: $ref_fasta"
        fi
    fi

    if [[ $errors -eq 0 ]]; then
        log_success "config.yaml validation passed"
        return 0
    else
        log_error "config.yaml validation failed with $errors error(s)"
        return 1
    fi
}

validate_sample_metadata() {
    local species=$1
    local species_dir="$BASE_DIR/$species"
    local metadata_file="$species_dir/sample_metadata.csv"
    local errors=0

    log_info "Validating sample metadata for: $species"

    # Check file exists
    if [[ ! -f "$metadata_file" ]]; then
        log_error "sample_metadata.csv not found: $metadata_file"
        return 1
    fi

    # Check file is not empty
    if [[ ! -s "$metadata_file" ]]; then
        log_error "sample_metadata.csv is empty"
        return 1
    fi

    # Check header format
    local header=$(head -n1 "$metadata_file")

    # Check if comma-separated (not tab-separated)
    if [[ "$header" =~ $'\t' ]]; then
        log_error "sample_metadata.csv appears to be tab-separated (should be comma-separated)"
        ((errors++))
    fi

    # Check required columns exist
    local required_cols=("sampID" "condition" "replicate")
    for col in "${required_cols[@]}"; do
        if ! echo "$header" | grep -q "$col"; then
            log_error "Missing required column: $col"
            ((errors++))
        fi
    done

    # Check age_category capitalization if present
    if echo "$header" | grep -q "age_category"; then
        local age_vals=$(tail -n +2 "$metadata_file" | cut -d, -f5 | sort -u)
        if echo "$age_vals" | grep -q "[a-z]"; then
            log_warning "age_category values should be uppercase (Y/O not y/o)"
        fi
    fi

    # Count samples
    local sample_count=$(tail -n +2 "$metadata_file" | grep -v '^[[:space:]]*$' | wc -l)
    log_info "Found $sample_count samples in metadata"

    if [[ $errors -eq 0 ]]; then
        log_success "sample_metadata.csv validation passed"
        return 0
    else
        log_error "sample_metadata.csv validation failed with $errors error(s)"
        return 1
    fi
}

count_samples() {
    local species=$1
    local species_dir="$BASE_DIR/$species"
    local samples_dir="$species_dir/samples"
    local count=0
    local valid_count=0

    if [[ ! -d "$samples_dir" ]]; then
        log_error "Samples directory not found: $samples_dir"
        return 0
    fi

    # Count sample directories
    while IFS= read -r -d '' sample_dir; do
        ((count++))
        # Check if FASTQ files exist
        if compgen -G "$sample_dir/*_1.fq.gz" > /dev/null && \
           compgen -G "$sample_dir/*_2.fq.gz" > /dev/null; then
            ((valid_count++))
        else
            log_warning "$(basename "$sample_dir"): No paired FASTQ files found"
        fi
    done < <(find "$samples_dir" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

    if [[ $count -eq 0 ]]; then
        log_error "No sample directories found in: $samples_dir"
        return 0
    fi

    log_info "Found $valid_count valid samples (out of $count directories)" >&2
    echo "$valid_count"
}

deploy_species() {
    local species=$1
    local cpus=${2:-$DEFAULT_CPUS}
    local mem_per_cpu=${3:-$DEFAULT_MEM_PER_CPU}
    local species_dir="$BASE_DIR/$species"

    log_section "Deploying Pipeline for $species"

    # Run validation first
    local validation_passed=true

    if ! validate_species_config "$species"; then
        validation_passed=false
    fi

    if ! validate_sample_metadata "$species"; then
        validation_passed=false
    fi

    if [[ "$validation_passed" == "false" ]]; then
        log_error "Validation failed. Fix errors before deployment."
        return 1
    fi

    # Count samples
    local sample_count=$(count_samples "$species")
    if [[ $sample_count -eq 0 ]]; then
        log_error "No valid samples found. Cannot deploy."
        return 1
    fi

    # Calculate resources
    local total_mem=$((cpus * mem_per_cpu))
    log_info "Array size: 1-${sample_count} (${sample_count} samples)"
    log_info "Resources per task: ${cpus} CPUs, ${mem_per_cpu}GB/CPU (${total_mem}GB total)"
    log_info "Estimated total: $((sample_count * cpus)) CPUs, $((sample_count * total_mem))GB"

    # Ask for confirmation
    echo ""
    echo -e "${YELLOW}Two-Job Workflow:${NC}"
    echo "  Job 1: Array job (process ${sample_count} samples in parallel)"
    echo "  Job 2: Analysis job (peak calling, deepTools, R analysis)"
    echo ""
    echo -e "${YELLOW}Ready to submit jobs. Press Enter to continue or Ctrl+C to cancel...${NC}"
    read -r

    # ===== JOB 1: Submit array job for sample processing =====
    log_info "Submitting array job for sample processing..."

    local array_job_cmd="sbatch --array=1-${sample_count} \
        --cpus-per-task=${cpus} \
        --mem-per-cpu=${mem_per_cpu}G \
        --export=SPECIES_NAME=${species},R_MODULE=${R_MODULE} \
        '$PIPELINE_SCRIPT'"

    local array_job_output=$(eval "$array_job_cmd")
    local array_job_id=$(echo "$array_job_output" | grep -oP 'Submitted batch job \K[0-9]+')

    if [[ -z "$array_job_id" ]]; then
        log_error "Array job submission failed"
        return 1
    fi

    log_success "Array job submitted: $array_job_id"

    # ===== JOB 2: Submit analysis job (depends on array completion) =====
    log_info "Submitting analysis job (depends on array completion)..."

    # Analysis job needs more memory for peak calling and R analysis
    local analysis_cpus=8
    local analysis_mem_per_cpu=8

    local analysis_job_cmd="sbatch --dependency=afterok:${array_job_id} \
        --cpus-per-task=${analysis_cpus} \
        --mem-per-cpu=${analysis_mem_per_cpu}G \
        --export=SPECIES_NAME=${species},R_MODULE=${R_MODULE} \
        '$PIPELINE_SCRIPT'"

    local analysis_job_output=$(eval "$analysis_job_cmd")
    local analysis_job_id=$(echo "$analysis_job_output" | grep -oP 'Submitted batch job \K[0-9]+')

    if [[ -z "$analysis_job_id" ]]; then
        log_error "Analysis job submission failed"
        log_warning "Array job $array_job_id is still running, but no follow-up job scheduled"
        return 1
    fi

    log_success "Analysis job submitted: $analysis_job_id"

    # ===== Display summary =====
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}Jobs submitted successfully!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    log_info "Species: $species"
    log_info "Job 1 (Array): $array_job_id (${sample_count} tasks)"
    log_info "Job 2 (Analysis): $analysis_job_id (depends on Job 1)"
    echo ""
    echo "Monitor jobs:"
    echo "  squeue -u \$USER"
    echo "  sacct -j ${array_job_id},${analysis_job_id} --format=JobID,JobName,State,ExitCode,Elapsed"
    echo ""
    echo "View array job logs:"
    echo "  tail -f atac-seq_pipeline-${array_job_id}_*.SLURMout"
    echo ""
    echo "View analysis job log:"
    echo "  tail -f atac-seq_pipeline-${analysis_job_id}.SLURMout"
    echo ""
    echo "Workflow:"
    echo "  1. Array tasks process samples (trim, align) → exit"
    echo "  2. Analysis job waits for all array tasks to complete"
    echo "  3. Analysis job runs peak calling, deepTools, R analysis"
    echo ""
    return 0
}

validate_all() {
    log_section "Validating All Species"

    local species_count=0
    local valid_count=0
    local invalid_count=0

    # Find all species directories
    for species_dir in "$BASE_DIR"/*/; do
        if [[ -f "$species_dir/config.yaml" ]]; then
            species=$(basename "$species_dir")
            ((species_count++))

            echo ""
            log_info "Checking: $species"

            if validate_species_config "$species" && validate_sample_metadata "$species"; then
                local sample_count=$(count_samples "$species")
                log_success "$species: Ready ($sample_count samples)"
                ((valid_count++))
            else
                log_error "$species: Failed validation"
                ((invalid_count++))
            fi
        fi
    done

    echo ""
    log_section "Validation Summary"
    echo "Total species: $species_count"
    echo "Valid: $valid_count"
    echo "Invalid: $invalid_count"

    return 0
}

deploy_all() {
    local cpus=${1:-$DEFAULT_CPUS}
    local mem_per_cpu=${2:-$DEFAULT_MEM_PER_CPU}

    log_section "Deploying All Species"

    # First, collect and validate all species
    local -a valid_species=()
    local -a valid_sample_counts=()
    local species_count=0

    echo ""
    log_info "Scanning for species in: $BASE_DIR"
    echo ""

    for species_dir in "$BASE_DIR"/*/; do
        if [[ -f "$species_dir/config.yaml" ]]; then
            local species=$(basename "$species_dir")
            ((species_count++))

            log_info "Validating: $species"

            if validate_species_config "$species" && validate_sample_metadata "$species"; then
                local sample_count=$(count_samples "$species")
                if [[ $sample_count -gt 0 ]]; then
                    valid_species+=("$species")
                    valid_sample_counts+=("$sample_count")
                    log_success "$species: Valid ($sample_count samples)"
                else
                    log_warning "$species: No valid samples found, skipping"
                fi
            else
                log_error "$species: Validation failed, skipping"
            fi
        fi
    done

    if [[ ${#valid_species[@]} -eq 0 ]]; then
        log_error "No valid species found to deploy"
        return 1
    fi

    # Show summary and confirm
    echo ""
    log_section "Deployment Summary"
    echo "Total species found: $species_count"
    echo "Valid species to deploy: ${#valid_species[@]}"
    echo ""
    echo "Species list:"
    for i in "${!valid_species[@]}"; do
        echo "  - ${valid_species[$i]} (${valid_sample_counts[$i]} samples)"
    done
    echo ""
    echo -e "${YELLOW}Ready to deploy ${#valid_species[@]} species. Press Enter to continue or Ctrl+C to cancel...${NC}"
    read -r

    # Deploy each species
    local -a all_job_ids=()
    local deployment_count=0

    for i in "${!valid_species[@]}"; do
        local species="${valid_species[$i]}"
        local sample_count="${valid_sample_counts[$i]}"

        echo ""
        log_section "Deploying: $species ($((i+1))/${#valid_species[@]})"

        # Submit array job
        local array_job_cmd="sbatch --array=1-${sample_count} \
            --cpus-per-task=${cpus} \
            --mem-per-cpu=${mem_per_cpu}G \
            --export=SPECIES_NAME=${species},R_MODULE=${R_MODULE} \
            '$PIPELINE_SCRIPT'"

        local array_job_output=$(eval "$array_job_cmd")
        local array_job_id=$(echo "$array_job_output" | grep -oP 'Submitted batch job \K[0-9]+')

        if [[ -z "$array_job_id" ]]; then
            log_error "Array job submission failed for $species"
            continue
        fi

        log_success "Array job submitted: $array_job_id"

        # Submit analysis job
        local analysis_cpus=8
        local analysis_mem_per_cpu=8

        local analysis_job_cmd="sbatch --dependency=afterok:${array_job_id} \
            --cpus-per-task=${analysis_cpus} \
            --mem-per-cpu=${analysis_mem_per_cpu}G \
            --export=SPECIES_NAME=${species},R_MODULE=${R_MODULE} \
            '$PIPELINE_SCRIPT'"

        local analysis_job_output=$(eval "$analysis_job_cmd")
        local analysis_job_id=$(echo "$analysis_job_output" | grep -oP 'Submitted batch job \K[0-9]+')

        if [[ -z "$analysis_job_id" ]]; then
            log_error "Analysis job submission failed for $species"
            log_warning "Array job $array_job_id is running without follow-up"
            all_job_ids+=("$species:$array_job_id:FAILED")
        else
            log_success "Analysis job submitted: $analysis_job_id"
            all_job_ids+=("$species:$array_job_id:$analysis_job_id")
            ((deployment_count++))
        fi
    done

    # Final summary
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}Deployment Complete${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    log_info "Successfully deployed: $deployment_count/${#valid_species[@]} species"
    echo ""
    echo "Job Summary:"
    for job_info in "${all_job_ids[@]}"; do
        IFS=':' read -r species array_id analysis_id <<< "$job_info"
        if [[ "$analysis_id" == "FAILED" ]]; then
            echo "  $species: Array=$array_id, Analysis=FAILED"
        else
            echo "  $species: Array=$array_id, Analysis=$analysis_id"
        fi
    done
    echo ""
    echo "Monitor all jobs:"
    echo "  squeue -u \$USER"
    echo ""
    echo "Check specific species logs:"
    for job_info in "${all_job_ids[@]}"; do
        IFS=':' read -r species array_id analysis_id <<< "$job_info"
        echo "  # $species"
        echo "  tail -f atac-seq_pipeline-${array_id}_*.SLURMout"
        if [[ "$analysis_id" != "FAILED" ]]; then
            echo "  tail -f atac-seq_pipeline-${analysis_id}.SLURMout"
        fi
    done
    echo ""

    return 0
}

# Main script logic
main() {
    local action=""
    local species=""
    local cpus=$DEFAULT_CPUS
    local mem=$DEFAULT_MEM_PER_CPU
    local validate_only=false

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --species)
                species="$2"
                action="deploy"
                shift 2
                ;;
            --deploy-all)
                action="deploy-all"
                shift
                ;;
            --validate)
                validate_only=true
                shift
                ;;
            --cpus)
                cpus="$2"
                shift 2
                ;;
            --mem)
                mem="$2"
                shift 2
                ;;
            --help)
                usage
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                usage
                exit 1
                ;;
        esac
    done

    # Check prerequisites
    if ! check_prerequisites; then
        exit 1
    fi

    # Execute action
    if [[ "$validate_only" == "true" ]]; then
        if [[ -n "$species" ]]; then
            # Validate single species
            if validate_species_config "$species" && validate_sample_metadata "$species"; then
                local sample_count=$(count_samples "$species")
                log_success "Validation passed for $species ($sample_count samples)"
                exit 0
            else
                log_error "Validation failed for $species"
                exit 1
            fi
        else
            # Validate all species
            validate_all
            exit 0
        fi
    elif [[ "$action" == "deploy" ]]; then
        if [[ -z "$species" ]]; then
            log_error "No species specified. Use --species NAME"
            usage
            exit 1
        fi
        deploy_species "$species" "$cpus" "$mem"
        exit $?
    elif [[ "$action" == "deploy-all" ]]; then
        deploy_all "$cpus" "$mem"
        exit $?
    else
        log_error "No action specified"
        usage
        exit 1
    fi
}

# Run main function
main "$@"

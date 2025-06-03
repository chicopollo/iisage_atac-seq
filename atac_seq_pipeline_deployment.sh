#!/bin/bash

####################################################################################
# Script: atac_seq_pipeline_deployment.sh
# Purpose: Helper script to deploy ATAC-seq pipeline jobs
# Author: Louis Paul Decena-Segarra, PhD
# Created: 2025-06-02
# Usage Examples:
#   ./deploy_atac_pipeline.sh --all                    # Deploy all species
#   ./deploy_atac_pipeline.sh --species chrysemys_picta # Deploy single species
#   ./deploy_atac_pipeline.sh --list                   # List available species
#   ./deploy_atac_pipeline.sh --validate               # Validate all configurations
####################################################################################

set -euo pipefail

# Configuration
BASE_DIR="/mnt/scratch/decenalo/atac_seq_projects"
PIPELINE_SCRIPT="/path/to/multi_species_atac_pipeline.sh"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --all                Deploy pipeline for all species"
    echo "  --species NAME       Deploy pipeline for specific species"
    echo "  --list               List all available species"
    echo "  --validate           Validate all species configurations"
    echo "  --status             Check status of running jobs"
    echo "  --setup SPECIES      Create directory structure for new species"
    echo "  --help               Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 --all"
    echo "  $0 --species chrysemys_picta"
    echo "  $0 --setup new_species_name"
}

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

check_prerequisites() {
    # Check if base directory exists
    if [[ ! -d "$BASE_DIR" ]]; then
        log_error "Base directory not found: $BASE_DIR"
        exit 1
    fi

    # Check if pipeline script exists
    if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
        log_error "Pipeline script not found: $PIPELINE_SCRIPT"
        exit 1
    fi

    # Check if SLURM is available
    if ! command -v sbatch >/dev/null 2>&1; then
        log_error "SLURM (sbatch) not available"

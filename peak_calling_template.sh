#!/bin/bash --login
##### SLURM Resource Requests #####

#SBATCH --time=48:00:00                          #How long the job will run (days-hours:minutes)
#SBATCH --nodes=1                               #How many compute nodes the job needs
#SBATCH --ntasks=1                              #How many concurrent tasks the job needs
#SBATCH --cpus-per-task=16                       #How many CPUs each task needs
#SBATCH --mem-per-cpu=8G                        #How much memory each CPU needs

##### SLURM Administrative Settings #####

#SBATCH --job-name 2025_ATACseq_cpicta_peak_calling                   #Name the job for convenience
#SBATCH --output=%x-%j.SLURMout                 #Name the output file (JobName-JobNumber.SLURMout)
#SBATCH --mail-type=ALL                         #Tell SLURM to email you when job starts, stops, error
#SBATCH --mail-user=decenalo@msu.edu                            #Provide SLURM your email address
#SBATCH --error=%x-%j.SLURMerr


####################################################################################
# Script: peak_calling.sh
# Purpose: Bulk ATAC-seq peak calling with size selection
# Author: Louis Paul Decena-Segarra, PhD
# Created: 2025-05-29
# Usage: sbatch peak_calling.sh
# Notes: Customize variables section before running!!
###################################################################################

set -euo pipefail  # Exit on error, undefined vars, pipe failures

##### bash  Commands to Run #####
#unload all modules
module purge

#Change directory to main working area in hpcc
cd /mnt/scratch/decenalo/

#Load required modules
#TROUBLESHOOTING
#If this doesn't work properly, make sure you are loading the correct modules
#use `$module spider GCC` (or any module) and check your versions are up to date and you are loading any pre-reqs
#module load GCC/13.2.0
#module load CUDA/12.4.0
#module load OpenMPI/4.1.6-GCC-13.2.0
module load MACS2/2.2.9.1-foss-2022b
module load SAMtools/1.17-GCC-12.2.0
# Define directories
PROCESSED_DIR=/mnt/scratch/decenalo/processed_files
PEAKS_DIR=/mnt/scratch/decenalo/peaks_output

# Create output directory if it doesn't exist
mkdir -p $PEAKS_DIR

# The GENOME_SIZE variable is used in macs2 callpeak function. This value can be found in the NCBI's genome assembly page
# Adjust based on the species you are working with
#Chrysemys_picta
GENOME_SIZE="2.4e9"

# ONCE ALL VARIABLES ARE DEFINED THERE SHOULDN'T BE ANY REASON TO MESS WITH THE FOLLOWING nucleotide
# THIS SHOULD WORK WITH ANY VOLUME OF SAMPLES, AS LONG AS EACH SAMPLE'S SORTED bam FILES ARE IN THE PROCESSED_DIR directory
# ADJUST THE NUMBER OF CPUs ACCORDING TO YOUR DATA VOLUME
# Ensure the script exits upon first error
set -e

# Loop through the BAM files in the processed directory

for BAM_FILE in $PROCESSED_DIR/*_sorted.bam
do
    SAMPLE_NAME=$(basename $BAM_FILE _sorted.bam)

    # Select reads <= 100 bp, exclude unwanted chromosomes, and create a new filtered BAM
    samtools view -@ 16 -h $BAM_FILE | \
    awk -v LEN=100 '{if ($9 <= LEN && $9 >= -LEN && $9 != 0 || $1 ~ /^@/) print $0}' | \
    grep -v "_random\|chrUn\|chrM" | \
    samtools view -@ 16 -uh -q 5 - | \
    samtools sort -@ 16 -o ${PEAKS_DIR}/${SAMPLE_NAME}_size_selected_100.bam
    samtools index -@ 16 ${PEAKS_DIR}/${SAMPLE_NAME}_size_selected_100.bam

    # Call peaks using MACS2
    macs2 callpeak -t ${PEAKS_DIR}/${SAMPLE_NAME}_size_selected_100.bam \
                   -f BAM \
                   -g $GENOME_SIZE \
                   --nomodel \
                   --nolambda \
                   --keep-dup all \
                   --call-summits \
                   -n $SAMPLE_NAME \
                   --outdir $PEAKS_DIR

    echo "Peaks called for $SAMPLE_NAME"
done

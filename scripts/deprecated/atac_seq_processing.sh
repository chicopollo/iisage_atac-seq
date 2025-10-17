#!/bin/bash --login

##### SLURM Resource Requests #####

#SBATCH --time=48:00:00                          #How long the job will run (days-hours:minutes)
#SBATCH --nodes=1                               #How many compute nodes the job needs
#SBATCH --ntasks=1                              #How many concurrent tasks the job needs
#SBATCH --cpus-per-task=16                       #How many CPUs each task needs
#SBATCH --mem-per-cpu=8G                        #How much memory each CPU needs

##### SLURM Administrative Settings #####

#SBATCH --job-name 2025_Cpicta_atac-seq                   #Name the job for convenience
#SBATCH --output=%x-%j.SLURMout                 #Name the output file (JobName-JobNumber.SLURMout)
#SBATCH --mail-type=ALL                         #Tell SLURM to email you when job starts, stops, error
#SBATCH --mail-user=decenalo@msu.edu                            #Provide SLURM your email address
#SBATCH --error=%x-%j.SLURMerr


####################################################################################
# Script: atac_seq_processing.sh
# Purpose: Bulk ATAC-seq trimming, alignment, and BAM generation
# Author: Louis Paul Decena-Segarra, PhD
# Created: 2025-05-29
# Usage: sbatch atac_seq_processing.sh
# Notes: Customize variables section before running!!
###################################################################################



##### bash  Commands to Run #####
#unload all modules
module purge

#Change directory to main working area in hpcc
cd /mnt/scratch/decenalo/
#Load required modules
#TROUBLESHOOTING
#If this doesn't work properly, make sure you are loading the correct modules
#use `$module spider Trimmomatic` (or any module) and check your versions are up to date and you are loading any pre-reqs
module load Trimmomatic
module load Bowtie2
module load samtools

################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Update these paths according to your setup!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
################################################

#The SAMPLES_DIR variable tells the script where are your fastq samples
SAMPLES_DIR=/mnt/scratch/decenalo/atac_seq
#The PROCESSED_DIR variable tells the script where to put your processed files
PROCESSED_DIR=/mnt/scratch/decenalo/processed_files
#The REF_GENOME varable tells the script where the genome index can be found (or created in case it doesn't exist)
REF_GENOME=/mnt/scratch/decenalo/genome_index_c_picta
#The REF_GENOME_FASTA tells the script where the fasta (nucleotide) sequence can be found
REF_GENOME_FASTA=/mnt/scratch/decenalo/cpicta_assembly/GCF_011386835.1_ASM1138683v2_genomic.fna
#The ADAPTERS_PATH variable points to the adapter list to be trimmed from fastq files
ADAPTERS_PATH=/mnt/home/decenalo/Documents/adapters/novogene.fa
#The TRIM_SETTINGS variable indicates the trimming settings for trimmomatic
TRIM_SETTINGS="ILLUMINACLIP:$ADAPTERS_PATH:2:30:10 SLIDINGWINDOW:4:15 MINLEN:25"


# ONCE ALL VARIABLES ARE DEFINED THERE SHOULDN'T BE ANY REASON TO MESS WITH THE NEXT PART OF THE SCRIPT
# THIS SHOULD WORK WITH ANY VOLUME OF SAMPLES, AS LONG AS EACH SAMPLE FASTQ FILES ARE IN THEIR OWN DIRECTORY

# Check if Bowtie2 index exists, if not build it
if [ ! -f "$REF_GENOME.1.bt2" ]; then
    echo "Index not found, creating Bowtie2 index..."
    bowtie2-build $REF_GENOME_FASTA $REF_GENOME
    echo "Index created."
else
    echo "Bowtie2 index found, proceeding..."
fi

# Process each sample directory
for SAMPLE_DIR in $SAMPLES_DIR/*; do
    SAMPLE_NAME=$(basename $SAMPLE_DIR)
    echo "Processing $SAMPLE_NAME"

    # Merge gzipped FASTQ files
    zcat $SAMPLE_DIR/*_1.fq.gz > $PROCESSED_DIR/${SAMPLE_NAME}_R1_merged.fastq
    zcat $SAMPLE_DIR/*_2.fq.gz > $PROCESSED_DIR/${SAMPLE_NAME}_R2_merged.fastq

    # Trim adapters and low quality
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 16 \
        $PROCESSED_DIR/${SAMPLE_NAME}_R1_merged.fastq \
        $PROCESSED_DIR/${SAMPLE_NAME}_R2_merged.fastq \
        $PROCESSED_DIR/${SAMPLE_NAME}_R1_paired.fastq $PROCESSED_DIR/${SAMPLE_NAME}_R1_unpaired.fastq \
        $PROCESSED_DIR/${SAMPLE_NAME}_R2_paired.fastq $PROCESSED_DIR/${SAMPLE_NAME}_R2_unpaired.fastq \
        $TRIM_SETTINGS

# Remove merged FASTQ files to save space (MSU's hpcc scratch autodeletes files after 45 days, so space is not really a concern right now)
# If you are processing a high volume of samples, uncomment the next two lines
#  rm $PROCESSED_DIR/${SAMPLE_NAME}_R1_merged.fastq
#  rm $PROCESSED_DIR/${SAMPLE_NAME}_R2_merged.fastq

    # Alignment with Bowtie2
    bowtie2 -x $REF_GENOME -1 $PROCESSED_DIR/${SAMPLE_NAME}_R1_paired.fastq \
            -2 $PROCESSED_DIR/${SAMPLE_NAME}_R2_paired.fastq \
            -S $PROCESSED_DIR/${SAMPLE_NAME}_aligned.sam \
            --very-sensitive --no-mixed --no-discordant --dovetail -X 2000 -k 2 -p 16

    # Convert SAM to BAM, sort and index
    samtools view -bS $PROCESSED_DIR/${SAMPLE_NAME}_aligned.sam | samtools sort -o $PROCESSED_DIR/${SAMPLE_NAME}_sorted.bam
    samtools index $PROCESSED_DIR/${SAMPLE_NAME}_sorted.bam

# Clean up unpaired files to save space (MSU's hpcc scratch autodeletes files after 45 days, so space is not really a concern right now)
# If you are processing a large quantity of data, uncomment the next 4 lines
  #  rm $PROCESSED_DIR/${SAMPLE_NAME}_R1_paired.fastq
  #  rm $PROCESSED_DIR/${SAMPLE_NAME}_R2_paired.fastq
  #  rm $PROCESSED_DIR/${SAMPLE_NAME}_R1_unpaired.fastq
  #  rm $PROCESSED_DIR/${SAMPLE_NAME}_R2_unpaired.fastq
done

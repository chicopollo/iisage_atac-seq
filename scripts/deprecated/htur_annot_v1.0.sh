#!/bin/bash --login
##### SLURM Resource Requests #####

#SBATCH --time=96:00:00                          #How long the job will run (days-hours:minutes)
#SBATCH --nodes=1                               #How many compute nodes the job needs
#SBATCH --ntasks=1                              #How many concurrent tasks the job needs
#SBATCH --cpus-per-task=32                       #How many CPUs each task needs
#SBATCH --mem-per-cpu=8G                        #How much memory each CPU needs

##### SLURM Administrative Settings #####

#SBATCH --job-name htur_annot                   #Name the job for convenience
#SBATCH --output=%x-%j.SLURMout                 #Name the output file (JobName-JobNumber.SLURMout)
#SBATCH --mail-type=ALL                         #Tell SLURM to email you when job starts, stops, error
#SBATCH --mail-user=decenalo@msu.edu                            #Provide SLURM your email address
#SBATCH --error=%x-%j.SLURMerr


####################################################################################
# Script: htur_annot.sh
# Purpose: quick and dirty genome annotation
# Author: Louis Paul Decena-Segarra, PhD
# Created: 2025-07-15
# Usage: sbatch htur_annot.sh
# Notes: Customize variables section before running!!
###################################################################################

module purge
module load Miniforge3            # or module load conda/…
conda activate funannote


# point Funannotate at your DB
export FUNANNOTATE_DB=/mnt/scratch/decenalo/annotation/funannotate_db
mkdir -p $FUNANNOTATE_DB


# Variables — adjust as needed
GENOME="H_turcicus_M_TG4449_HiFiasmHiC.no_mtDNA_p_ctg.FINAL_SL.fasta"
SPECIES="Hemidactylus turcicus"            # for metadata (no quotes in funannotate call)
BUSCO_DB="vertebrata_odb10"             # change to your lineage
SLURM_CPUS_PER_TASk=32
CPUS=$SLURM_CPUS_PER_TASK
prefix=$(basename "$GENOME")

# Make all output dirs
mkdir -p clean_out sorted_out masked_out predict_out annotate_out


# 1) Clean (remove redundant scaffolds) — output is a FASTA file
funannotate clean \
  -i "$GENOME" \
  -o "clean_out/$prefix" \
  --cpus "$CPUS"

# 2) Sort (by length & rename headers) — output is a FASTA file
funannotate sort \
  -i "clean_out/$prefix" \
  -o "sorted_out/$prefix"

# 3) Mask repeats (RepeatModeler + RepeatMasker) — output dir
funannotate mask \
  -i "sorted_out/$prefix" \
  -o "masked_out/$prefix" \
  --cpus $CPUS

# 4) Predict genes (ab initio + BUSCO training) — output dir
funannotate predict \
  -i "masked_out/$prefix" \
  -o "predict_out/$prefix" \
  --species "$SPECIES" \
  --busco_db "$BUSCO_DB" \
  --cpus "$CPUS"

# 5) Functional annotation — output dir
funannotate annotate \
  -i predict_out \
  -o annotate_out \
  --cpus "$CPUS"

conda deactivate

# IISAGE: Multi-Species ATAC-seq pipeline
Comprehensive and standarized ATAC-seq analysis and figure generation pipeline for SLURM based HPCC environments.

## Quick start
1. Clone repository
2. Organize data in the described directory structure
3. Edit config/config.yaml (1 file for each species) and place it in the species folder
4. Adjust paths and module names to match your environment in pipeline scripts
5. Submit job using SBATCH ('multi_species_atac_pipeline.sh')
6. Generate figures and tables using R ('generate_atac_figs.R')

Note: I will integrate the figure generation to the main pipeline so it will all come out automatically.

## Features
- Process multiple species at a time, ideal for comparative analysis
- Automated and standarized: QC, Trimming, alignment, peak calling, visualizations
- R based downstream analysis
- Modular stages (can be enabled/disabled independently)

## Requirements
Verify the following modules are available in your set up  by using `$module spider MODULE_NAME`. Check closely what are the versions available, how to properly load them, and modify line 286 in the atac_seq_pipeline_deployment.sh script.

Trimmomatic
Bowtie2
SAMtools
MACS2
deepTools
R

## Directory structure

For this pipeline to be completed effectively, organize your files in the following manner.
```
 $BASE_DIR/
 ├── species1/
 │   ├── config.yaml
 │   ├── samples/
 │   │   ├── sample1/
 |   |       ├──Reads1_fwd.fastq.gz
 |   |       └──Reads2_rev_fastq.gz
 │   │   └── sample2/
 │   ├── reference/
 │   └── sample_metadata.xlsx
 └── species2/
     ├── config.yaml
     ├── samples/
     ├── reference/
     └── sample_metadata.xlsx
```
## Test dataset
In order to test this pipeline, you should send a message to abroniko@msu.edu to request a dataset of three reptile species (_Chrysemys picta_, _Thamnophis elegans_, and _Hemydactylus turcicus_). This dataset is already organized and contains the config.yaml files needed to run.

If you are part of the Bronikowski Lab, you can access the data through the following path in "Carson":
```
/data/grpdata/broniko_lab/IISAGE_atac_seq/test_data
```
Within this folder you will find subfolders for each species.

From your HPCC environment navigate to the directory in which you have the pipeline scripts and use rsync or scp to copy the data:

```
rsync -r $USER@carson.kbs.msu.edu:/data/grpdata/broniko_lab/IISAGE_atac_seq/test_data .
```
or
```
scp -r $USER@carson.kbs.msu.edu:/data/grpdata/broniko_lab/IISAGE_atac_seq/test_data .
```

## Troubleshooting
There are instances in which the modules of your particular HPCC environment might not work so well together, or in which additional modules might be requires. While editing the script for submission, it is a good practice to run `$module spider MODULE_NAME` to analyze what are the main requirements

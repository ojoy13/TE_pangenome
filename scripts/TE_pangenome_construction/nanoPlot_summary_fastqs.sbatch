#!/bin/bash

## Script for running NanoPlot on all fastq.gz files
## Note: NanoStat requires python 3
## Example usage:
## sbatch nanoPlot_summary.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J nanoPlot
#SBATCH -o /Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%A_%a.err

pwd; hostname; date

# Set query files
queries=($(ls $inDir/*fastq.gz| xargs -n 1 basename))

# Set variables
inDir="/TE_pangenome/TE_pangenome_construction/outputs/fastqs"
outDir="/TE_pangenome/TE_pangenome_construction/outputs/NanoPlot_fastqs"

# Make outDir, if it doesn't exist
mkdir -p $outDir

# Run NanoPlot
echo "Starting NanoPlot..."
echo "Processing all fastq.gz files in $inDir..."

# Run NanoPlot on all fastq.gz files
NanoPlot -t 2 --fastq $inDir/*.fastq.gz \
--maxlength 40000 \
--plots dot \
--legacy hex \
-o $outDir

echo $(date +"[%b %d %H:%M:%S] Done!")

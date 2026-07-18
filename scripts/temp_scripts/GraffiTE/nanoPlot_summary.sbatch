#!/bin/bash

## Script for running NanoPlot on all fastq.gz
##
## Example usage:
## inDir=. outDir=. sbatch nanoPlot_summary.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J nanoPlot
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%A_%a.err

# Activate virtual environment that contains NanoPlot
source ~/jhub_venv_3.10.8/bin/activate

# Set query files
queries=($(ls $inDir/*fastq.gz| xargs -n 1 basename))

# Make output dir, if it doesn't exist
mkdir -p $outDir

# Run NanoPlot
pwd; hostname; date

echo "Starting NanoPlot..."
echo "Processing all fastq.gz files in $inDir..."

# Run NanoPlot on all fastq.gz files
NanoPlot -t 2 --fastq $inDir/*.fastq.gz \
--maxlength 40000 \
--plots dot \
--legacy hex \
-o $outDir/combined_plots

echo $(date +"[%b %d %H:%M:%S] Done!")

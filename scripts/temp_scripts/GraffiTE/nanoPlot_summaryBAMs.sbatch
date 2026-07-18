#!/bin/bash

## Script for running NanoPlot on all .bam
##
## Example usage:
## inDir=. outDir=. sbatch nanoPlot_summaryBAMs.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J nanoPlot_BAM
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%A_%a.err

# Activate virtual environment that contains NanoPlot
source ~/jhub_venv_3.10.8/bin/activate

# Set query files
queries=($(ls $inDir/*.bam | xargs -n 1 basename))

# Make output dir, if it doesn't exist
mkdir -p $outDir

# Run NanoPlot
pwd; hostname; date

echo "Starting NanoPlot..."
echo "Processing all bam files in $inDir..."

# Run NanoPlot on all .bam files
NanoPlot -t 12 \
--color green \
--bam $inDir/*.bam \
--downsample 10000 \
-o $outDir/bamplots_downsampled

echo $(date +"[%b %d %H:%M:%S] Done!")

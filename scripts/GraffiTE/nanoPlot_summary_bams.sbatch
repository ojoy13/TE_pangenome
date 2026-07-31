#!/bin/bash

## Script for running NanoPlot on all .bam
## Note: NanoStat requires python 3
## Example usage:
## sbatch nanoPlot_summaryBAMs.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J nanoPlot_BAM
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

# Set Variables
inDir="/TE_pangenome/GraffiTE/outputs/NanoStat_bams" 
outDir="/TE_pangenome/GraffiTE/outputs/NanoPlot_bams" 

# Make outDir, if it doesn't exist
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
-o $outDir

echo $(date +"[%b %d %H:%M:%S] Done!")

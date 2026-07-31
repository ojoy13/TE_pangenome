#!/bin/bash

## Script for running NanoStat on fastq.gz files
## Note: NanoStat requires python 3
## Example usage:
## sbatch nanostat.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --time=24:00:00
#SBATCH --mem=128GB
#SBATCH --array 0-147 

# Job name and output
#SBATCH -J nanoStat
#SBATCH -o /Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%A_%a.err

pwd; hostname; date

# Set query files
queries=($(ls $inDir/*.fastq.gz | xargs -n 1 basename))

# Set variables
inDir="/TE_pangenome/GraffiTE/outputs/fastqs"
outDir="/TE_pangenome/GraffiTE/outputs/NanoStat_fastqs"

# Make outDir, if it doesn't exist
mkdir -p $outDir

# Run hisat2
echo "[%b %d %H:%M:%S] Starting NanoStat..."
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

NanoStat --fastq $inDir/${queries[$SLURM_ARRAY_TASK_ID]} \
--name ${queries[$SLURM_ARRAY_TASK_ID]%.fastq}_stats.txt \
--outdir $outDir \
--threads 32

echo $(date +"[%b %d %H:%M:%S] Done!")
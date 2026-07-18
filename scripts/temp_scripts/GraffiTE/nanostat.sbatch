#!/bin/bash

## Script for running NanoStat
##
## Example usage:
## inDir=. outDir=nanostat_output sbatch --array 0-40 nanostat.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J nanoStat
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%A_%a.err

# Activate virtual environment that contains NanoStat
source ~/jhub_venv_3.10.8/bin/activate

# Set query files
queries=($(ls $inDir/*.fastq.gz | xargs -n 1 basename))

# Make output dir, if it doesn't exist
mkdir -p $outDir

# Run hisat2
pwd; hostname; date

echo "[%b %d %H:%M:%S] Starting NanoStat..."
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

NanoStat --fastq $inDir/${queries[$SLURM_ARRAY_TASK_ID]} \
--name ${queries[$SLURM_ARRAY_TASK_ID]%.fastq}_stats.txt \
--outdir $outDir \
--threads 32

echo $(date +"[%b %d %H:%M:%S] Done!")
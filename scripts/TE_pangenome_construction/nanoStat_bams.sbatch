#!/bin/bash

## Script for running NanoStat for merged and sorted bams
## Note: NanoStat requires python 3
## Example usage:
## sbatch nanostat_bam.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=50GB
#SBATCH --array 0-20

# Job name and output
#SBATCH -J nanoStat_bam
#SBATCH -o /Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%A_%a.err

pwd; hostname; date

# Set query files
queries=($(ls $inDir/*.bam | xargs -n 1 basename))


# Set Variables
inDir="/TE_pangenome/TE_pangenome_construction/outputs/bams_merged_sorted"
outDir="/TE_pangenome/TE_pangenome_construction/outputs/NanoStat_bams"

# Make outDir, if it doesn't exist
mkdir -p $outDir

# Run NanoStat
echo $(date + "[%b %d %H:%M:%S] Starting NanoStat...")
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

NanoStat --bam $inDir/${queries[$SLURM_ARRAY_TASK_ID]} \
--outdir $outDir \
--name ${queries[$SLURM_ARRAY_TASK_ID]%.merged.bam}.txt
--threads 8

echo $(date +"[%b %d %H:%M:%S] Done!")

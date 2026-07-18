#!/bin/bash

## Script for indexing bams
##
## Example usage:
## inDir=. outDir=. sbatch --array 0-20 index_bams.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=50GB

# Job name and output
#SBATCH -J index_bams
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%A_%a.err

module load samtools/1.8

#set query files
queries=($(ls $inDir/*.bam | xargs -n 1 basename))

echo "[%b %d %H:%M:%S] Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

samtools index -@ 8 -b $inDir/${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done!")
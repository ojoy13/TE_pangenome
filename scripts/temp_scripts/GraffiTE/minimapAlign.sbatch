#!/bin/bash

## Script for aligning long-read WGS bams
##
## Example usage:
## inDir=. outDir=. sbatch --array 0-147 minimapAlign.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=96:00:00
#SBATCH --mem=240GB

# Job name and output
#SBATCH -J minimapAlign
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%A_%a.err

#load modules
module load samtools
module load minimap2

#set query files
queries=($(ls $inDir/*.fastq.gz | xargs -n 1 basename))

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Aligning reads and sorting file...")

#start
minimap2 -t 16 -ax map-ont /Shares/CL_Shared/db/genomes/hg38/fa/hg38.main.fa $inDir/${queries[$SLURM_ARRAY_TASK_ID]} \
| samtools sort - -@16 -o $inDir/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Indexing bam file...")

samtools index -@16 $inDir/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Done!")

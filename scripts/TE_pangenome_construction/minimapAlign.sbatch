#!/bin/bash

## Script for aligning long-read WGS fastq.gz files from HGSVC3
##
## Example usage:
##  sbatch minimapAlign.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=96:00:00
#SBATCH --mem=240GB
#SBATCH --array 0-147 

# Job name and output
#SBATCH -J minimapAlign
#SBATCH -o /Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%A_%a.err

pwd; hostname; date

#set query files
queries=($(ls $inDir/*.fastq.gz | xargs -n 1 basename))

# set variables
inDir="/TE_pangenome/TE_pangenome_construction/outputs/fastqs"
outDir="/TE_pangenome/TE_pangenome_construction/outputs/bams_bySRA"
ref="/path/to/dir/hg38.main" 

# Make outDir, if it doesn't exist
mkdir -p $outDir

# run Minimap2 align
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Aligning reads and sorting file...")

#start
minimap2 -t 16 -ax map-ont $ref $inDir/${queries[$SLURM_ARRAY_TASK_ID]} \
| samtools sort - -@16 -o $inDir/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Indexing bam file...")

samtools index -@16 $inDir/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Done!")

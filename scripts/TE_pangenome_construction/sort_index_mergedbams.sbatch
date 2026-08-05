#!/bin/bash

## Script for running merge_bams.py to merge bams into samples IDs
##
## Example usage:
## sbatch  merge_bams.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=72:00:00
#SBATCH --mem=50GB
#SBATCH --array=0-20

# Job name and output
#SBATCH -J sort_index_merge_bams
#SBATCH -o /Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%A_%a.err

pwd; hostname; date

# Set variables
inDir="/TE_pangenome/TE_pangenome_construction/outputs/bams_merged"
outDir="/TE_pangenome/TE_pangenome_construction/outputs/bams_merged_sorted"

# Make outDir, if it doesn't exist
mkdir -p $outDir

#set query files
queries=($(ls $inDir/*.merged.bam | xargs -n 1 basename))

# run sort & index merged bam
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Aligning reads and sorting file...")

samtools sort -@16 -o $outDir/${queries[$SLURM_ARRAY_TASK_ID]%.merged.bam}.merged.sorted.bam \
$inDir/${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Indexing bam file...")

samtools index -@16 $outDir/${queries[$SLURM_ARRAY_TASK_ID]%.merged.bam}.merged.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Done!")

echo "Done processing sample: ${queries[$SLURM_ARRAY_TASK_ID]%.merged.bam}"
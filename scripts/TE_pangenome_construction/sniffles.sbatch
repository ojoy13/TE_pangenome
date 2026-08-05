#!/bin/bash

## Script for identifying structural variants by SAMPLE
## Note: Sniffles2 requires python 3
## Example usage:
## sbatch sniffles.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=4:00:00
#SBATCH --mem=10GB
#SBATCH --array 0-20

# Job name and output
#SBATCH -J sniffles
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

pwd; hostname; date

#set query files
queries=($(ls $inDir/*.sorted.bam | xargs -n 1 basename))

# Set variables
inDir="/TE_pangenome/TE_pangenome_construction/outputs/bams_merged_sorted"
outDir="/TE_pangenome/TE_pangenome_construction/outputs/sniffles"
ref="/path/to/dir/hg38.main" 

# Make outDir, if it doesn't exist
mkdir -p $outDir

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Starting sniffles...")

#start
sniffles --minsvlen 100 --threads 16 --reference $ref \
--input $inDir/${queries[$SLURM_ARRAY_TASK_ID]} --snf $outDir/${queries[$SLURM_ARRAY_TASK_ID]%.merged.sorted.bam}.snf \
--vcf $outDir/${queries[$SLURM_ARRAY_TASK_ID]%.merged.sorted.bam}.vcf

echo $(date +"[%b %d %H:%M:%S] Done!")
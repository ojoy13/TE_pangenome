#!/bin/bash

## Script for identifying structural variants by SAMPLE
##
## Example usage:
## inDir=bams_withID outDir=bams_withID sbatch --array 0-20 sniffles.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=4:00:00
#SBATCH --mem=10GB

# Job name and output
#SBATCH -J sniffles
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%A_%a.err

#activate virtual environment that contains sniffles
source ~/jhub_venv_3.10.8/bin/activate

#set query files
queries=($(ls $inDir/*.bam | xargs -n 1 basename))

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Starting sniffles...")

#start
sniffles --minsvlen 100 --threads 16 --reference /Shares/CL_Shared/db/genomes/hg38/fa/hg38.main.fa --input $inDir/${queries[$SLURM_ARRAY_TASK_ID]} --snf $inDir/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.snf --vcf $inDir/${queries[$SLURM_ARRAY_TASK_ID]%.sorted.bam}.vcf

echo $(date +"[%b %d %H:%M:%S] Done!")
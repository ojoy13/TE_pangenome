#!/bin/bash

## Script to gzip fastqs for quality control report
##
## Example usage:
## sbatch --array 0-45 gzip.sbatch

# General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --time=12:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J gzip
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%A_%a.err


#set variables:
inDir="/scratch/Users/olde5615/data/pangenome.v3_RNA/HGSVC_shortread_RNA_fastqs"

echo "Gzipping file: "${queries[$SLURM_ARRAY_TASK_ID]}
queries=($(ls $inDir/*.fastq | xargs -n 1 basename))

gzip -v -f ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done!")

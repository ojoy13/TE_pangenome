#!/bin/bash

## Script for compressing HGSVC3 fastq SRA files into gzip format
## Example usage:
## sbatch gzip.sbatch

# General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=72:00:00
#SBATCH --mem=10GB
#SBATCH --array=0-147

# Job name and output
#SBATCH -J gzip_fastq
#SBATCH -o /Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%A_%a.err

pwd; hostname; date

#set variables
inDir="/TE_pangenome/TE_pangenome_construction/outputs/fastqs"

#set query files
queries=($(ls $inDir/*.fastq | xargs -n 1 basename))

echo $(date +"[%b %d %H:%M:%S] Gzip compressing ${queries[$SLURM_ARRAY_TASK_ID]}...")

# run
gzip $inDir/${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done!")

#!/bin/bash

## Example usage:
## sbatch --array=0-40 gzip.sbatch

# General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=72:00:00
#SBATCH --mem=10GB

# Job name and output
#SBATCH -J gzip_fastq
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%A_%a.err

inDir="/scratch/Users/olde5615/data/scale_up"

#set query files
queries=($(ls $inDir/*.fastq | xargs -n 1 basename))

echo $(date +"[%b %d %H:%M:%S] Gzip compressing ${queries[$SLURM_ARRAY_TASK_ID]}...")

# run
gzip $inDir/${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done!")

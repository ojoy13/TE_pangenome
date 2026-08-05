#!/bin/bash

## Script for running deeptools bamCoverage
##
## Example usage:
## sbatch bam_to_bw.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=8:00:00
#SBATCH --mem=32GB
#SBATCH --array 0-21

# Job name and output
#SBATCH -J bam_to_bw
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

pwd; hostname; date

# Set variables
inDir="/TE_pangenome/outputs/RNA/hisat2"  
outDir="/TE_pangenome/outputs/RNA/bigWigs" 
binSize=1
numCPU=8

# Make outDir, if it doesn't exist
mkdir -p $outDir

# Set query files
queries=($(ls $inDir/*.uniq.bam | xargs -n 1 basename))

echo "Starting bam2bigwig forward..."

# Forward bigwig
singularity exec --bind /path/to/dir/singularity/deeptools-3.0.1-py35_1.img bamCoverage \
-b $inDir/${queries[$SLURM_ARRAY_TASK_ID]} \
-o $outDir/${queries[$SLURM_ARRAY_TASK_ID]%.uniq.bam}_fwd.bw \
--binSize=$binSize -p $numCPU --normalizeUsing CPM \
--filterRNAstrand forward --ignoreForNormalization chrX chrM

echo "Starting bam2bigwig reverse..."
# Reverse bigwig
singularity exec --bind /path/to/dir/singularity/deeptools-3.0.1-py35_1.img bamCoverage \
-b $inDir/${queries[$SLURM_ARRAY_TASK_ID]} \
-o $outDir/${queries[$SLURM_ARRAY_TASK_ID]%.uniq.bam}_rev.bw \
--binSize=$binSize -p $numCPU --normalizeUsing CPM \
--filterRNAstrand reverse --ignoreForNormalization chrX chrM

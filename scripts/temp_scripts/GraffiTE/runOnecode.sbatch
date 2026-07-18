#!/bin/bash

## Script for running OneCodeToFindThemAll
##
## Example usage: 
## sbatch runOnecode.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J runOnecode
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%j.err


FASTA_FILE=/scratch/Users/olde5615/data/scale_up/bams_withID/repeatmasker_dir/indels.fa
REPMASK_DIR=/scratch/Users/olde5615/data/scale_up/bams_withID/repeatmasker_dir
REPMASK_OUT=${REPMASK_DIR}/$(basename ${FASTA_FILE}).out

### Step 2 run OneCodeToFindThemAll
echo "Running onecode..."
one_code_to_find_them_all.pl --rm ${REPMASK_OUT} --unknown --ltr ${REPMASK_DIR}/OneCode_LTR.dic > ${REPMASK_DIR}/onecode.log 2>&1

echo $(date +"[%b %d %H:%M:%S] Done!")

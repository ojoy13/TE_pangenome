#!/bin/bash

## Script for builing onecode LTR library
##
## Example usage:
## sbatch buildOnecode.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J buildOnecode
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%j.err


FASTA_FILE=/scratch/Users/olde5615/data/scale_up/bams_withID/repeatmasker_dir/indels.fa
REPMASK_DIR=/scratch/Users/olde5615/data/scale_up/bams_withID/repeatmasker_dir
REPMASK_OUT=${REPMASK_DIR}/$(basename ${FASTA_FILE}).out

### Step 1 make the dictionary file if LTRs are split between _LTR and _I
echo "Building onecode LTR dictionary..."
build_dictionary.pl --rm ${REPMASK_OUT} --unknown > ${REPMASK_DIR}/OneCode_LTR.dic

echo "Done building onecode LTR dictionary"
echo $(date +"[%b %d %H:%M:%S] Done!")

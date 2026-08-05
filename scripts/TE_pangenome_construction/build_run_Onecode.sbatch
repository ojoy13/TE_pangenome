#!/bin/bash

## Script for building and running Onecode LTR library
##
## Example usage:
## sbatch buildOnecode.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J build_run_Onecode
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

pwd; hostname; date

# Set variables
REPMASK_DIR="/TE_pangenome/TE_pangenome_construction/outputs/repeatmasker"
FASTA_FILE="$REPMASK_DIR/indels.fa"
REPMASK_OUT=${REPMASK_DIR}/$(basename ${FASTA_FILE}).out

# Make the dictionary file if LTRs are split between _LTR and _I
echo "Building onecode LTR dictionary..."
build_dictionary.pl --rm ${REPMASK_OUT} --unknown > ${REPMASK_DIR}/OneCode_LTR.dic

echo "Done building onecode LTR dictionary"
echo $(date +"[%b %d %H:%M:%S] Done!")

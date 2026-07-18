#!/bin/bash

## Script for Filtering Repeats for LINE-1 & SVA-VNTR
##
##
## Example usage:
## sbatch filterRepeats.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J filterRepeats.sbatch
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%j.err


FASTA_FILE=/scratch/Users/olde5615/data/scale_up/bams_withID/repeatmasker_dir/indels.fa
REPMASK_DIR=/scratch/Users/olde5615/data/scale_up/bams_withID/repeatmasker_dir
REPMASK_OUT=${REPMASK_DIR}/$(basename ${FASTA_FILE}).out

REPMASK_ONECODE_OUT=${REPMASK_DIR}/$(basename ${FASTA_FILE}).onecode.out

echo $(date +"[%b %d %H:%M:%S] Filtering onecode repeats started!")


### Step 3 parse the output
# remove output if previous run
echo "Concatenate outputs..."
rm ${REPMASK_DIR}/ALL.elem_sorted.csv 2>/dev/null

# loop over each individual output and concatenate in a single file (can't cat, often there are too many files)
for i in ${REPMASK_DIR}/*.elem_sorted.csv
do
    cat $i >> ${REPMASK_DIR}/ALL.elem_sorted.csv
done

# parse the concatenated outputs and save the raw to recover stitched copies
echo "Parse outputs..."
cat <(head -n 3 ${REPMASK_OUT} | sed 's/ID/ID\tfragmts/g') <(grep '###' ${REPMASK_DIR}/ALL.elem_sorted.csv | cut -f 1-16 | sed 's/###//g') > ${REPMASK_ONECODE_OUT}
mv ${REPMASK_DIR}/ALL.elem_sorted.csv ${REPMASK_DIR}/ALL.onecode.elem_sorted.bak

# clean-up the mess (BE SURE YOU HAVE RENAMED THE OUTPUT)
echo "Cleanup..."
find ${REPMASK_DIR} -name "*.copynumber.csv" -delete
find ${REPMASK_DIR} -name "*.elem_sorted.csv" -delete
find ${REPMASK_DIR} -name "*.ltr.csv" -delete
find ${REPMASK_DIR} -name "*.transposons.csv" -delete

echo $(date +"[%b %d %H:%M:%S] Filtering Done!")

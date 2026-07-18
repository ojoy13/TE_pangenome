#!/bin/bash

## Script for Repeat Masking a VCF file
##
## Example usage for vcf containing 5 individuals increase time for bigger population:
## inFile=genotypes.vcf outFile=genotypes_repmasked.vcf.gz sbatch repmask_vcf.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=100:00:00
#SBATCH --mem=500GB

# Job name and output
#SBATCH -J repmask_vcf
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%j.err

#load modules
module load bcftools
module load samtools

# find all insertions and deletions and put into indels.fa
FASTA_FILE=indels.fa
bcftools view --types indels --include 'ILEN>0' ${inFile} | grep -v "#" | awk '{print(sprintf(">%s\n%s", $3, $5))}' > ${FASTA_FILE}
bcftools view --types indels --include 'ILEN<0' ${inFile} | grep -v "#" | awk '{print(sprintf(">%s\n%s", $3, $4))}' >> ${FASTA_FILE}

# verify that indels.fa ids are not longer than 50 characters
if grep ">" indels.fa | awk 'length > 50' | grep -q .; then
    echo "variant IDs must be no greater than 50 characters"
    exit 1
fi

mkdir repeatmasker_dir
REPMASK_DIR=repeatmasker_dir


# run RepeatMasker
echo "processing repeats"

RepeatMasker -nolow -lib /Shares/CL_Shared/db/repbase25.09/Libraries/RMRBSeqs.fasta -s -dir ${REPMASK_DIR} -pa 16 ${FASTA_FILE}

echo "RepeatMasker Done!"


REPMASK_OUT=${REPMASK_DIR}/$(basename ${FASTA_FILE}).out

echo $(date +"[%b %d %H:%M:%S] Done!")

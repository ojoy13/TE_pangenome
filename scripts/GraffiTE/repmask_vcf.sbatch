#!/bin/bash

## Script for Repeat Masking a VCF file
## Note: Download latest version - this script uses RepeatMasker v4.1.0
## Example usage for vcf containing 21 individuals:
## sbatch repmask_vcf.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=100:00:00
#SBATCH --mem=500GB

# Job name and output
#SBATCH -J repmask_vcf
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

pwd; hostname; date

# Set variables
inFile="/TE_pangenome/GraffiTE/outputs/sniffles_popn/sniffles2_variants.vcf" 
outDir="/TE_pangenome/GraffiTE/outputs/repeatmasker"

# Make outDir, if it doesn't exist
mkdir -p $outDir

# find all insertions and deletions and put into indels.fa
FASTA_FILE=indels.fa
bcftools view --types indels --include 'ILEN>0' ${inFile} | grep -v "#" | awk '{print(sprintf(">%s\n%s", $3, $5))}' > ${FASTA_FILE}
bcftools view --types indels --include 'ILEN<0' ${inFile} | grep -v "#" | awk '{print(sprintf(">%s\n%s", $3, $4))}' >> ${FASTA_FILE}

# verify that indels.fa ids are not longer than 50 characters
if grep ">" indels.fa | awk 'length > 50' | grep -q .; then
    echo "variant IDs must be no greater than 50 characters"
    exit 1
fi

# run RepeatMasker
echo "processing repeats"

RepeatMasker -nolow -s -species human -dir $outDir -pa 16 ${FASTA_FILE}

echo "RepeatMasker Done!"


REPMASK_OUT=$outDir/$(basename ${FASTA_FILE}).out

echo $(date +"[%b %d %H:%M:%S] Done!")

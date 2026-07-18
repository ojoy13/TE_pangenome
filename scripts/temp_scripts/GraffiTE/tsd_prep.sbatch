#!/bin/bash

## Script for preparing target site duplication identification
##
##
## Example usage:
## sbatch tsd_prep.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=50GB

# Job name and output
#SBATCH -J tsd_prep
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%j.err

# Load modules
module load bedtools
module load bcftools


# variable list
VCF=/scratch/Users/olde5615/data/scale_up/bams_withID/repeatmasker_dir/genotypes_repmasked_filtered.vcf # filtered vcf with repeatmasker
REF=/Shares/CL_Shared/db/genomes/hg38/fa/hg38.main.fa # ref genome
WIN=30 # windows size in flanking to search TSD - 30 is the GraffiTE default
MSK=/scratch/Users/olde5615/data/scale_up/bams_withID/repeatmasker_dir/indels.fa.masked # masked SVs from repeatmasker
OUT_VCF=/scratch/Users/olde5615/data/scale_up/bams_withID/TSD_output/pangenie.vcf
OUT_DIR=/scratch/Users/olde5615/data/scale_up/bams_withID/TSD_output

###################################################
# Step 1: extract flanking of each retained TE SV #
###################################################
echo "extracting flanking..."

mkdir -p ${OUT_DIR}

# get contig length for bedtools
bcftools view ${VCF} | grep "contig=" | sed 's/\#\#contig=<ID=//g;s/,length=/\t/g;s/>//g' > ${OUT_DIR}/gLength.txt

# create a bed with vcf entries 
bcftools view ${VCF} | grep -v '#' | \
 grep 'n_hits=1;\|n_hits=2;' | \
 grep -v 'mam_filter_2=VNTR_ONLY' | \
 awk '/n_hits=1/ && length($4) < length($5) {print $1"\t"$2"\t"($2)+1"\t"$3; next} /n_hits=1/ && length($4) > length($5) {print $1"\t"$2"\t"($2+length($4))"\t"$3; next} /n_hits=2/ && length($4) < length($5) && /5P_INV/ {print $1"\t"$2"\t"($2)+1"\t"$3; next} /n_hits=2/ && length($4) > length($5) && /5P_INV/ {print $1"\t"$2"\t"($2+length($4))"\t"$3}' > ${OUT_DIR}/oneHit_SV_coordinates.bed

echo $(date +"[%b %d %H:%M:%S] Created bedfile with vcf entries")


# extend +/- ${WIN} bp in two entries per SV
echo $(date +"[%b %d %H:%M:%S] Starting SV flanking +/-30bp")

cat <(bedtools slop -i ${OUT_DIR}/oneHit_SV_coordinates.bed -g ${OUT_DIR}/gLength.txt -l 30 -r 0 | awk '{print $0"__L"}') \
<(bedtools slop -i ${OUT_DIR}/oneHit_SV_coordinates.bed -g ${OUT_DIR}/gLength.txt -l 0 -r 30 | awk '{print $0"__R"}') | \
sort -k1,1 -k2,2n -k3,3n | awk '/__L/ {print $1"\t"$2"\t"($2+30)"\t"$4; next} /__R/ {print $1"\t"($3-30)"\t"$3"\t"$4}' > ${OUT_DIR}/oneHit_SV_coordinates_win.bed
# extract fasta from flanking
bedtools getfasta -fi ${REF} -bed ${OUT_DIR}/oneHit_SV_coordinates_win.bed -name > ${OUT_DIR}/flanking_sequences.fasta

echo $(date +"[%b %d %H:%M:%S] Done flanking")


##################################################
# Step 2: extract 5' and 3' of each masked TE SV #
##################################################
echo "extracting SVs' 5' and 3' ends..."

# filter the indels.fa.masked to keep only single RM hits (1 TE per SV)
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(cut -f 4 ${OUT_DIR}/oneHit_SV_coordinates.bed) ${MSK} > ${OUT_DIR}/oneHit_indels.fa.masked
# linearize fasta, then trim and split in two seq (L and R)
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' ${OUT_DIR}/oneHit_indels.fa.masked | \
awk -v len=${WIN} -F '\t' '{x=len;L=length($2);printf("%s\n%s\n%s\n%s\n",$1"__L",(L<=x?$2:substr($2,2,x+1)),$1"__R",(L<=x?$2:substr($2,1+L-x,x)));}' > ${OUT_DIR}/SV_sequences_L_R_trimmed_WIN.fa
# export the list of SV to search TSD for next process parallelization
grep '>' ${OUT_DIR}/SV_sequences_L_R_trimmed_WIN.fa | sed 's/>//g;s/__/\t/g' | cut -f 1 | sort | uniq > ${OUT_DIR}/indels.txt

echo $(date +"[%b %d %H:%M:%S] Done!")

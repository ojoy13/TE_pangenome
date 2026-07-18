#!/bin/bash

## Script for annotating repeat files
##
##
## Example usage:
## sbatch annotateRepeats.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J annotateRepeats
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%j.err

# Load modules
module load samtools/1.8
module load bedtools
module load bcftools

FASTA_FILE=/scratch/Users/olde5615/data/scale_up/bams_withID/repeatmasker_dir/indels.fa
REPMASK_DIR=/scratch/Users/olde5615/data/scale_up/bams_withID/repeatmasker_dir
REPMASK_OUT=${REPMASK_DIR}/$(basename ${FASTA_FILE}).out
REPMASK_ONECODE_OUT=${REPMASK_DIR}/$(basename ${FASTA_FILE}).onecode.out
VCF=/scratch/Users/olde5615/data/scale_up/bams_withID/sniffles2_variants.vcf

#part 4 annotation vcf
#create empty file
ANNOT_FILE=${REPMASK_DIR}/vcf_annotation

# vcf_annotation_1 is a temporary created from annotate_vcf.R

# calculate the total span of TEs on the SV without overlap
echo $(date +"[%b %d %H:%M:%S] Starting...")
echo "compute repeat proportion for each SVs..."
samtools faidx ${REPMASK_DIR}/indels.fa

awk '{print $1"\t"$2}' indels.fa.fai > ${REPMASK_DIR}/indels.length
awk 'NR > 3 {print $5"\t"$6"\t"$7"\t"$10}' ${REPMASK_ONECODE_OUT} | bedtools merge > ${REPMASK_DIR}/merge.bed

RMQUERIES=$(awk 'NR > 3 {print $5}' ${REPMASK_ONECODE_OUT} | sort | uniq)
rm -rf span &> /dev/null # clean in case there is a "span" file already
rm -rf ${ANNOT_FILE}.gz &> /dev/null # clean in case there was a ${ANNOT_FILE}.gz file already

##replace with
for i in ${RMQUERIES}; do
    # Extract matching line from merge.bed
    merge_result=$(grep -w "${i}" ${REPMASK_DIR}/merge.bed | awk '{print ($3-$2)}' | paste -sd+ | bc)
    
    # Extract matching line from indels.length
    indel_result=$(grep -w "${i}" ${REPMASK_DIR}/indels.length)

    # Check if both merge_result and indel_result are not empty
    if [[ -z "$merge_result" || -z "$indel_result" ]]; then
        echo -e "${i}\tNA\tNA\tNA" >> ${REPMASK_DIR}/span
    else
        # Extract the fourth column from indel_result (length)
        length=$(echo "$indel_result" | awk '{print $2}')
        
        # Check if length is zero to avoid division by zero error
        if [[ "$length" -eq 0 ]]; then
            echo -e "${i}\t${merge_result}\t${length}\tNA" >> ${REPMASK_DIR}/span
        else
            ratio=$(echo "scale=6; ${merge_result} / ${length}" | bc)
            echo -e "${i}\t${merge_result}\t${length}\t${ratio}" >> ${REPMASK_DIR}/span
        fi
    fi
done

#get lines without NA - span 18055 lines & span_no_NA 11031 lines
grep -v 'NA' ${REPMASK_DIR}/span > ${REPMASK_DIR}/span_no_NA

# merge with ${ANNOT_FILE}
join -13 -21 -a1 <(sort -k3,3 ${ANNOT_FILE}_1)  <(sort -k1,1 ${REPMASK_DIR}/span_no_NA) | sed 's/ /\t/g' | \
 awk '{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15}' | \
 awk '{if (NF == 10) {print $0"\t0\t0"} else {print $0}}' | \
 sort -k1,1 -k2,2n > ${ANNOT_FILE}

# For mammal filtered files
 echo "Mammalian filters ON. Filtering..."
 
 # FILTER 1: L1 Twin Priming and similar
 TwP=$(awk '{if ($6 == 2 && $11 == "C,+" && $10~/LINE/) {names = split($10,n,",",seps); ids = split($12, i,",",seps); if (names n[1] == names n[2] && ids i[1] == ids i[2]) {print $0"\tTwP"}}}' ${ANNOT_FILE} | awk '{print $3}')
 rm ${REPMASK_DIR}/TwP.txt &> /dev/null
 
 for i in $TwP 
 do
 grep -w $(grep -w "${i}" <(grep -v "#" ${VCF} | awk '{print $1"_"$2"\t"$3}') | cut -f 2) ${REPMASK_DIR}/indels.fa.onecode.out | \
 awk 'getline second {line=split(second, a, "\t", sep); print $12"\t"$13"\t"$14"\t"a[12]"\t"a[13]"\t"a[14]}' | \
 awk -v coord=${i} '{if ($2<$5) {print coord"\t5P_INV:plus"} else {print coord"\t5P_INV:minus"}}' >> ${REPMASK_DIR}/TwP.txt
 done

 # FILTER 2: SVA VNTR
 # get coordinates of each putative TE in TE (single TE hits)
 awk '$6 == 1 && $14 > 0.9 {print $1"_"$2"\t"$1"\t"$2"\t"($2+1)"\t"$8"\t"$10}' ${ANNOT_FILE} | grep 'SVA' | sort -k1,1 > ${REPMASK_DIR}/SVA_candidates
 join -17 -25 <(join -11 -21 SVA_candidates <(awk 'NR > 3 && !/^#/ {print $1"_"$2"\t"$3}' ${VCF} | \
 sort -k1,1) | sort -k7,7) <(sort -k5,5 ${REPMASK_DIR}/indels.fa.onecode.out) | \
 sed 's/ /\t/g' | \
 awk '{if ($15 == "+") {print $16"\t"$18"\t"$19"\t"$1} else {print $16"\t"$20"\t"$19"\t"$1}}' > ${REPMASK_DIR}/SVA_candidates.bed
 
 # create bed with SVAs' VNTR position according to TRF on DFAM3.6 models (2022)
 rm ${REPMASK_DIR}/SVA_VNTR.bed &> /dev/null

 echo -e "SVA_A\t436\t855\t+" >> ${REPMASK_DIR}/SVA_VNTR.bed
 echo -e "SVA_B\t431\t867\t+" >> ${REPMASK_DIR}/SVA_VNTR.bed
 echo -e "SVA_C\t432\t851\t+" >> ${REPMASK_DIR}/SVA_VNTR.bed
 echo -e "SVA_D\t432\t689\t+" >> ${REPMASK_DIR}/SVA_VNTR.bed
 echo -e "SVA_E\t428\t864\t+" >> ${REPMASK_DIR}/SVA_VNTR.bed
 echo -e "SVA_F\t435\t857\t+" >> ${REPMASK_DIR}/SVA_VNTR.bed
 bedtools intersect -wao -a <(sort -k1,1 -k2,2n ${REPMASK_DIR}/SVA_candidates.bed) -b SVA_VNTR.bed | awk '{if ($9/($3-$2) > 0.9) {print $4"\tVNTR_ONLY:"$5":"$2":"$3}}' | sort -k1,1 > ${REPMASK_DIR}/SVA_VNTR.txt

 join -a1 -13 -21 <(sort -k3,3  ${ANNOT_FILE}) <(sort -k1,1 human_pangenome/repeatmasker_dir/TwP.txt) | sed 's/ /\t/g' | awk '{if (NF == 14) {print $0"\tNone"} else {print $0}}' > ${REPMASK_DIR}/vcf_annotation.temp
 join -a1 -11 -21 <(sort -k1,1 ${REPMASK_DIR}/vcf_annotation.temp) ${REPMASK_DIR}/SVA_VNTR.txt | sed 's/ /\t/g' | awk '{if (NF == 15) {print $0"\tNone"} else {print $0}}' | \
 awk '{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' | sort -k1,1 -k2,2n > ${ANNOT_FILE}

 echo "writing vcf..."
 bgzip ${ANNOT_FILE}
 tabix -s1 -b2 -e2 ${ANNOT_FILE}.gz

 HDR_FILE=${REPMASK_DIR}/hdr_file

 echo -e '##INFO=<ID=n_hits,Number=1,Type=Integer,Description="Number of repeats found in insertion">' >> ${HDR_FILE}
 echo -e '##INFO=<ID=match_lengths,Number=.,Type=Integer,Description="Insertion lengths spanned by each repeat">' >> ${HDR_FILE}
 echo -e '##INFO=<ID=repeat_ids,Number=.,Type=String,Description="Repeat family IDs">' >> ${HDR_FILE}
 echo -e '##INFO=<ID=matching_classes,Number=.,Type=String,Description="Repeat class names">' >> ${HDR_FILE}
 echo -e '##INFO=<ID=fragmts,Number=.,Type=Integer,Description="Number of fragments merged into one by one code">' >> ${HDR_FILE}
 echo -e '##INFO=<ID=RM_hit_strands,Number=.,Type=String,Description="RepeatMasker hit strands">' >> ${HDR_FILE}
 echo -e '##INFO=<ID=RM_hit_IDs,Number=.,Type=String,Description="RepeatMasker hit IDs">' >> ${HDR_FILE}
 echo -e '##INFO=<ID=total_match_length,Number=1,Type=Integer,Description="Insertion length spanned by repeats">' >> ${HDR_FILE}
 echo -e '##INFO=<ID=total_match_span,Number=1,Type=Float,Description="Insertion span spanned by repeats">' >> ${HDR_FILE}
 echo -e '##INFO=<ID=mam_filter_1,Number=1,Type=String,Description="L1 with 5P inversion; strand">' >> ${HDR_FILE}
 echo -e '##INFO=<ID=mam_filter_2,Number=1,Type=String,Description=">90% of the SV is a SVA VNTR; SVA_family:consensus_hit_coordinates">' >> ${HDR_FILE}
 echo -e '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> ${HDR_FILE}

 cat <(grep '#' ${VCF}) <(grep -v '#' ${VCF} | sort -k1,1 -k2,2n) > ${REPMASK_DIR}/genotypes.sorted.vcf
 bcftools annotate -a ${ANNOT_FILE}.gz -h ${HDR_FILE} \
 -c CHROM,POS,ID,REF,ALT,INFO/n_hits,INFO/fragmts,INFO/match_lengths,INFO/repeat_ids,INFO/matching_classes,INFO/RM_hit_strands,INFO/RM_hit_IDs,INFO/total_match_length,INFO/total_match_span,INFO/mam_filter_1,INFO/mam_filter_2 genotypes.sorted.vcf | \
 bcftools view -Oz -o ${REPMASK_DIR}/genotypes_repmasked.vcf.gz

echo $(date +"[%b %d %H:%M:%S] Done!")

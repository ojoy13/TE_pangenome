#!/bin/bash

## Script for final report of target site duplications
##
##
## Example usage:
## sbatch tsd_report.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J tsd_report
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

pwd; hostname; date

# Set variables 
OUT_DIR="/TE_pangenome/outputs/TE_pangenome_construction/filtered_pangenome/TSD"
VCF="/TE_pangenome/outputs/TE_pangenome_construction/filtered_pangenome/pangenome3TE_TSDfilt_prep.vcf"
TSD_FILE="$OUT_DIR/TSD_annotation"
HDR_FILE="$OUT_DIR/header_file"
tsd_summary="$OUT_DIR/indels.txt.TSD_summary.txt"
tsd_full_log="$OUT_DIR/indels.txt.TSD_full_log.txt"
OUT_VCF="/TE_pangenome/outputs/TE_pangenome_construction/filtered_pangenome/TE_pangenome_21samples.vcf"

#rename tsd files
cat ${tsd_summary} > TSD_summary.txt
cat ${tsd_full_log} > TSD_full_log.txt


echo $(date +"[%b %d %H:%M:%S] Starting report...")

# combine summary with filtered repmasked genotypes
# corrected combine summary
join -13 -21 <(grep -v "#" ${VCF} | cut -f 1-5 | sort -k3,3) <(grep 'PASS' ${TSD_summary} | awk '{print $1"\t"$(NF-2)","$(NF-1)}' | sort -k1,1) | awk '{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$6}' | sort -k1,1 -k2,2n > ${TSD_FILE}
  
#add header line
echo -e '##INFO=<ID=TSD,Number=1,Type=String,Description="Target site duplication sequence passing filters">' >> ${HDR_FILE}
  
# zip and index
bgzip ${TSD_FILE}
tabix -f -s1 -b2 -e2 ${TSD_FILE}.gz
bcftools annotate -a ${TSD_FILE}.gz -h ${HDR_FILE} -c CHROM,POS,REF,ALT,INFO/TSD ${VCF} | bcftools view > $OUT_VCF
bgzip $OUT_VCF
tabix -p vcf ${OUT_DIR}.gz
  
echo $(date +"[%b %d %H:%M:%S] Done!")

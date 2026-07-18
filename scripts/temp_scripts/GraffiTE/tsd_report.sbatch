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
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%j.err

# Load modules
module load bcftools
module load samtools/1.8

# variable list
VCF=/scratch/Users/olde5615/data/scale_up/bams_withID/repeatmasker_dir/genotypes_repmasked_filtered.vcf
OUT_DIR=/scratch/Users/olde5615/data/scale_up/bams_withID/TSD_output

tsd_summary=/scratch/Users/olde5615/data/scale_up/bams_withID/TSD_output/indels.txt.TSD_summary.txt
tsd_full_log=/scratch/Users/olde5615/data/scale_up/bams_withID/TSD_output/indels.txt.TSD_full_log.txt

#rename tsd files
cat ${tsd_summary} > TSD_summary.txt
cat ${tsd_full_log} > TSD_full_log.txt


echo $(date +"[%b %d %H:%M:%S] Starting report...")

TSD_FILE=/scratch/Users/olde5615/data/scale_up/bams_withID/TSD_output/TSD_annotation


#combine summary with filtered repmasked genotypes
# corrected combine summary
join -13 -21 <(grep -v "#" ${VCF} | cut -f 1-5 | sort -k3,3) <(grep 'PASS' ${TSD_summary} | awk '{print $1"\t"$(NF-2)","$(NF-1)}' | sort -k1,1) | awk '{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$6}' | sort -k1,1 -k2,2n > ${TSD_FILE}

 
HDR_FILE=header_file
  
#add header line
echo -e '##INFO=<ID=TSD,Number=1,Type=String,Description="Target site duplication sequence passing filters">' >> ${HDR_FILE}
  
  
# zip and index
bgzip ${TSD_FILE}
tabix -f -s1 -b2 -e2 ${TSD_FILE}.gz
bcftools annotate -a ${TSD_FILE}.gz -h ${HDR_FILE} -c CHROM,POS,REF,ALT,INFO/TSD ${VCF} | bcftools view > TE_pangenome_21samples.vcf
  
echo $(date +"[%b %d %H:%M:%S] Done!")

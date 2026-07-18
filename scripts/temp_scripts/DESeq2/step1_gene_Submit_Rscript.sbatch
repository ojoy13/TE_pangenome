#!/bin/bash
#SBATCH --job-name=snif_transcript
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=olde5615@fiji.colorado.edu                     # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=04:00:00                                 # Time limit hrs:min:sec
#SBATCH --partition=short                               # Partition/queue requested on server    
#SBATCH --mem=3GB                      
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%j.err

################################################################                                                                                                                      
##Load required modules 
################################################################
module load R/4.4.0

################################################################ 
##Set paths
################################################################
rscript_path=/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts/deseq_1sniffle_forlooping.R

################################################################ 
##Job information
################################################################ 
printf "\nR code being run: ${rscript_path}"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou have requested $SLURM_CPUS_ON_NODE core(s).\n"

################################################################                  
## Run Script 
################################################################
#using the Rscript command to run the script
Rscript ${rscript_path} $a_sniffle 

echo Job finished at `date +"%T %a %d %b %Y"`

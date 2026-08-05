import os

# Run the script
command = "python3 fix_vcf.py --ref /TE_pangenome/GraffiTE/data/hg38.main.fa --vcf_in /TE_pangenome/GraffiTE/outputs/filtered_pangenome/pangenome3TE.vcf --vcf_out /TE_pangenome/GraffiTE/outputs/filtered_pangenome/pangenome3TE_TSDfilt.vcf"

exit_code = os.system(command)

if exit_code == 0:
    print("fix_vcf.py completed successfully!")
else:
    print(f"fix_vcf.py failed with exit code: {exit_code}")
import pandas as pd
import subprocess

# Not all samples were genotyped for each chromosome
# Create a list of fully genotyped samples

# Input and output files
INPUT_PATH="/path/to/input/files/merge" # Use the merged VCFs output from script 1_run_mergeSTR_statSTR.sh
OUTPUT_PATH="/path/to/output/sample/lists/"
OUTPUT_LIST="/path/to/final/output/sample/list.txt"

# Get samples genotyped for each chromosome
chromsomes = [str(s) for s in range(1, 23)] + ['X']

for chrom in chromsomes:
	command = f'bcftools query -l ${INPUT_PATH}/chr{chrom}.vcf.gz | sort | uniq > ${OUTPUT_PATH}/chr{chrom}.txt'
	subprocess.run(command, shell=True)

# Get samples genotypes for all chromosomes
full_list=[]
for chrom in chromsomes:
	filename = f'${OUTPUT_PATH}/chr{chrom}.txt'
	with open(filename, 'r') as f:
		list_for_chrom = f.readlines()
	list_for_chrom = [s.strip() for s in list_for_chrom]
	if full_list=[]:
		full_list=list_for_chrom
	else:
		full_list = list(set(full_list) & set(list_for_chrom))

# Save
fout=open(OUTPUT_LIST, 'w')
fout.write('\n'.join(full_list))
fout.close()

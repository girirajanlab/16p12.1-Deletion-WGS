import sys
import csv

# Parse PennCNV output

# Input and output files
indiv_calls="/path/to/individual/CNV/calls.txt" # Use the "penncnv_all.txt" output from script 1_run_PennCNV.sh here
family_calls="/path/to/family/CNV/calls.txt" # Use the "penncnv_family.txt" output from script 1_run_PennCNV.sh here
OUTPUT="/paht/to/output/file.txt"

# Open files
indiv=open(indiv_calls, 'rU')
family=open(family_calls, 'rU')
outwriter = csv.writer(file(OUTPUT,'w'),delimiter='\t')

# Read calls
all_lines=indiv.readlines()
family_lines=family.readlines()

# Iterate through individual CNV calls to merge individual and family calls together
merge_list=[]

for all_line in all_lines: #Iterate through individual lines
	dup_found=False
	all_line=all_line.split()
	for family_line in family_lines:
		family_line_split=family_line.split()
		family_test=family_line_split[:7] # Check if individual call matches family call
		if all_line == family_test: # Duplicate calls in individual and family records
			merge_list.append(family_line_split)
			dup_found=True
			family_lines.remove(family_line)
	if dup_found==False:
		merge_list.append(all_line) # No duplicate line found

for family_line in family_lines: # Add family-only lines
	merge_list.append(family_line.split())

# Parse each record and export to file
for line in merge_list:
	chrom=line[0].split(":")[0]
	start=line[0].split(":")[1].split("-")[0]
	end=line[0].split(":")[1].split("-")[1]
	numsnp=line[1].strip("numsnp=")
	length=line[2].strip('"length=').replace(",","")
	state=line[3].split(",")[0]
	cnv_type="dup" # Duplication
	if int(line[3].split(",")[1].strip("cn=").strip('"'))<2: # Deletion
		cnv_type="del"
	sample=line[4]
	startsnp="rs"+line[5].strip("startsnp=")
	endsnp=line[6].strip("endsnp=")
	try: # Add family info
		family=line[7]
	except(IndexError): # For indivudal-only calls
		family=""
	try: # Same as above...
		inh=line[8]
	except(IndexError):
		inh=""

	outwriter.writerow([chrom]+[start]+[end]+[numsnp]+[length]+[state]+[cnv_type]+[sample]+[startsnp]+[endsnp]+[family]+[inh])

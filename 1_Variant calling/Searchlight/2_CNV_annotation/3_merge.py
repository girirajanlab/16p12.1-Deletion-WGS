import csv

# Merge adjacent PennCNV calls if they:
# 1. Overlap (any %)
# 2. Have a gap of less than 50kb AND <20% of the combined CNV length

# Then filter out CNVs less than 50kb in length or than encompass less than 5 SNPs

# This step will need to be run for each batch

# Input and output files
INPUT="/path/to/input/file.txt" # Use the output of script 2_parse_penncnv.py here
SNP_LIST="/path/to/snp/list.tsv" # Use the same list of SNP positions from script 1_run_PennCNV.sh
OUTPUT="/path/to/output/file.txt"

gap_limit=0.2
gap_length=50000

size_limit=50000
snp_limit=5

# Helper functions
def mergeOverlap(start1,end1,start2,end2,recip):
	overlap=0
	if end2>=end1: #Partial overlap
		overlap=end1-start2
	else: #Full overlap
		overlap=end2-start2
	if (end1-start1)*recip<=overlap and (end2-start2)*recip<=overlap: #Check to see if CNVs are reciprocal by defined threshold
		return True
	return False

def mergeAdjacent(start1,end1,start2,end2,gap_limit,gap_length):
	gap=start2-end1
	sum_size=(end1-start1)+(end2-start2)
	if gap <= gap_length and float(gap)/sum_size<=gap_limit:
		return True
	return False 

def lineTransform(line,start,end,startsnp,endsnp,chrom,snp_lines):
	chrom=chrom.strip("chr")
	line[1]=start
	line[2]=end
	line[3]=snpcalc(chrom,start,end,snp_lines)
	line[4]=end-start
	line[8]=startsnp
	line[9]=endsnp
	return line

def snpcalc(chrom,start,end,snp_lines):
	counter=0
	for record in snp_lines:
		record=record.split("\t")
		if record[1]==chrom and int(record[2])>=start and int(record[2])<=end:
			counter=counter+1
	return counter

# Open files
infile=open(INPUT,'rU')
snp_file=open(SNP_LIST,'rU')
outwriter = csv.writer(file(OUTPUT,'w'),delimiter='\t')

# Read through lines and combine CNVs
snp_lines=snp_file.readlines()
inlines=infile.readlines()

checklist=[]
sample_name=""

for line1str in inlines:
	line1=line1str.strip().split("\t")

	if line1[7]!=sample_name:
		sample_name=line1[7]
		print "Processing sample "+sample_name

	newline=line1
	start1=int(newline[1])
	end1=int(newline[2])
	id_tag1=line1[7]+line1[6]+line1[0]+line1[1]+line1[2]
	if id_tag1 in checklist:
		continue
	else:
		checklist.append(id_tag1)

	# From PennCNV log file, remove samples with low QC

	# Merge adjacent and overlapping CNVs
	for line2str in inlines:
	 	line2=line2str.strip().split("\t")
	 	start2=int(line2[1])
	 	end2=int(line2[2])
	 	id_tag2=line2[7]+line2[6]+line2[0]+line2[1]+line2[2]
		
	 	if line2[7]==line1[7] and line2[0]==line1[0] and line2[6]==line1[6]: #Get all lines with same sample, chromosome and type
			if id_tag2 not in checklist:

				if start1<=start2: #Line 1 is first

					if end1>=start2: #CNVs overlap, Line 1 is first
	 					if mergeOverlap(start1,end1,start2,end2,recip)==True: #Passes reciprocal
	 						if end2>=end1: #Change line for partial overlap
	 							newline=lineTransform(newline,start1,end2,newline[8],line2[9],line2[0],snp_lines)
	 						if len(line2)>len(newline): #Add family info if not present
	 							newline.append(line2[10])
	 							newline.append(line2[11])
	 						checklist.append(id_tag2)

	 				else: #Check for merging CNVs
	 					if mergeAdjacent(start1,end1,start2,end2,gap_limit,gap_length)==True: 
	 						newline=lineTransform(newline,start1,end2,newline[8],line2[9],line2[0],snp_lines)
	 						if len(line2)>len(newline): #Add family info if not present
	 							newline.append(line2[10])
	 							newline.append(line2[11])
	 						checklist.append(id_tag2)

	 			else: #Line 2 is first

	 				if end2>=start1: #CNVs overlap, Line 2 is first
	 					if mergeOverlap(start2,end2,start1,end1,recip)==True: #Passes reciprocal
	 						if end1>=end2: #Change line for partial overlap
	 							newline=lineTransform(newline,start2,end1,line2[8],newline[9],line2[0 ],snp_lines)
	 						else:
	 							newline=lineTransform(newline,start2,end2,line2[8],line2[9],line2[0],snp_lines)
	 						if len(line2)>len(newline): #Add family info if not present
	 							newline.append(line2[10])
	 							newline.append(line2[11])
	 						checklist.append(id_tag2)

					else: #Check for merging CNVs
	 					if mergeAdjacent(start2,end2,start1,end1,gap_limit,gap_length)==True:
	 						newline=lineTransform(newline,start2,end1,line2[8],newline[9],line2[0],snp_lines)
	 						if len(line2)>len(newline): #Add family info if not present
	 							newline.append(line2[10])
	 							newline.append(line2[11])
	 						checklist.append(id_tag2)

	#Filter CNVs out by size and number of SNPs
	if int(newline[4])>=size_limit and int(newline[3])>=snp_limit:
		outwriter.writerow(newline)

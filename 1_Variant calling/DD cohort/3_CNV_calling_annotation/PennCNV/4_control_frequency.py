import csv

# Input and output files
INPUT="/path/to/input/file.txt" # Use the output of script 3_merge.py here
OUTPUT="/path/to/output/file.txt"

# Reference files
CONTROL_CNVS="/path/to/frequency/data/from/control/microarray/cohort/calls.txt
# These data are from Coe et al. Nat Genet 2014 (https://pubmed.ncbi.nlm.nih.gov/25217958/) and contain microarray CNV calls from a control cohort

# Helper functions
def intervalsize(start,end):
	if int(start) > int(end):
		return int(start)-int(end)
	return int(end) - int(start)

def getOverlap(start1,end1,start2,end2): #returns reciprocol overlap if overlap, else returns 0
	recip = 0
	if int(start1)>int(end1):
		holder=start1
		start1=end1
		end1=holder
	if int(start2)>int(end2):
		holder=start2
		start2=end2
		end2=holder
	if int(start1) >= int(start2) and int(end1) >= int(end2) and int(start1)<=int(end2):
		recip = int(end2)-int(start1)
	if int(start1) <= int(start2) and int(end1) <= int(end2) and int(start2)<=int(end1):
		recip = int(end1)-int(start2)
	if int(start1) <= int(start2) and int(end1) >= int(end2):
		recip = int(end2)-int(start2)
	if int(start1) >= int(start2) and int(end1) <= int(end2):
		recip = int(end1)-int(start1)
	return recip

def reciprocal_overlap(start1,end1,start2,end2,recip):
	len1=intervalsize(start1,end1)
	len2=intervalsize(start2,end2)
	over_len=getOverlap(start1,end1,start2,end2)
	if len1*recip<=over_len and len2*recip<=over_len:
		return True
	return False

# Open files
infile=open(INPUT,'rU')
control_file=open(CONTROL_CNVS,'rU')
outwriter = csv.writer(file(OUTPUT,'w'), delimiter='\t')

# Add header to output
header=["Chromsome","Start","End","SNPs","Size","CNV_type","Patient_ID","StartSNP","EndSNP","Control_freq","Control_freq_label"]
outwriter.writerow(header)

# Iterate through calls and identify their frequency in a control cohort
batch_sample_id=""
for inline in inlines:
	inline=inline.strip().split("\t")
	chrom=inline[0]
	start=int(inline[1])
	end=int(inline[2])
	cnv_type=inline[5]
	if cnv_type=="dup":
		cnv_type=3
	elif cnv_type=="del":
		cnv_type=1

	if inline[6]!=batch_sample_id:
		batch_sample_id=inline[6]
		print "Processing sample "+batch_sample_id

	sample_list=[]

	for control_line in control_lines:
		control_line=control_line.strip().split("\t")
		control_chrom=control_line[0]
		control_start=int(control_line[1])
		control_end=int(control_line[2])
		control_type=int(control_line[4])
		sample_id=control_line[3]

		if cnv_type==control_type and chrom==control_chrom:
			overlap=reciprocal_overlap(start,end,control_start,control_end,recip)
			if overlap==True:
				sample_list.append(sample_id)

	list(set(sample_list)) #Get rid of duplicate samples
	control_count=len(sample_list)
	frequency=float(control_count)/18572 # 18572 unique individuals in file
	outline=inline[0:9]+[control_count]
	if frequency<0.001: 
		outline.append("rare")
	else:
		outline.append("common")
	outwriter.writerow(outline)

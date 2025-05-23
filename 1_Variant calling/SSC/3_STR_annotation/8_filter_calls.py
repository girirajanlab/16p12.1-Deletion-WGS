

# Input and output files
INPUT_TABLE="/path/to/input/table.tsv" # Use the output of script 7_combine_chromosomes.sh here
OUTPUT_TABLE="/path/to/output/table.tsv"

# Load files
fin=open(INPUT_TABLE, 'r')
fout=open(OUTPUT_TABLE, 'w')

# Iterate through lines and check for exonic variants
for line in fin:
	line = line
	
	# header
	if line.startswith('chrom'):
		fout.write(line)

	sline  = line[:-1].split('\t')
	sample = sline[2]
	func   = sline[13]
	gene   = sline[14]
	loeuf  = sline[-1]
	
	if func != 'exonic':
		continue
	
	fout.write(line)
fin.close()
fout.close()

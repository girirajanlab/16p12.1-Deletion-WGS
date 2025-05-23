import pandas as pd
import sys

# Manta reports multiple SV types - extract deletions and duplications
input_file='/path/to/input.vcf' # Use the VCF generated from 1_run_manta.sh here
output_file='/path/to/output_table.txt'

# Get index of header lines to skip
exclude = [i for i, line in enumerate(open(file, 'r')) if line.startswith('#')]

# Read file as dataframe
vcf = pd.read_csv(sys.argv[1], sep = '\t', skiprows = exclude, names = ['Chr','Pos', 'ID', 'Ref', 'Alt', 'Qual', 'FILTER', 'INFO', 'FORMAT', 'FORMAT_VALUES'])

# Make a new column for Type for easy filtering
vcf['Type'] = vcf.INFO.str.split('SVTYPE=', expand = True)[1].str.split(';', expand = True)[0]

# Remove calls that are not deletions or duplications
vcf[(vcf.Type=='DEL') | (vcf.Type=='DUP')].to_csv(output_file, sep = '\t', index = False)


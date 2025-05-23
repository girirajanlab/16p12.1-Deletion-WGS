import pandas as pd

# Filter PennCNV calls based on frequency

# Input and output files
INPUT="/path/to/input/file.csv" # Use the output from script 7_gnomad_frequency.py
PATHOGENIC="/paht/to/pathogenic/cnv.bed" # Use tht output from script 8_pathogenic_overlap.sh
OUTPUT="/path/to/output/file.txt"

# Load files
calls=pd.read_csv(INPUT)
patho_cnvs=pd.read_csv(PATHOGENIC, sep='\t',
					names = ['Chromosome', 'Start', 'End', 'Type', 'PatientID',
							'Pathogenic_Chr', 'Pathogenic_Start', 'Pathogenic_End', 'Pathogenic_Name']))

# Merge calls with pathogenic CNV annotations
calls=pd.merge(calls, patho_cnvs[['Chromosome', 'Start', 'End', 'Type', 'PatientID','NPathogenic_Name']],
				on = ['Chromosome', 'Start', 'End', 'Type', 'PatientID'], how = 'left')

# Replace empty Pathogenic names with .
calls.loc[calls.Pathogenic_Name.isnull(), "Pathogenic_Name"] = '.'

# Filter for intracohort frequency
calls = calls[calls.intracohort_count <= 10]

# gnomAD frequency filter (do not apply this to Pathogenic CNVs)
calls = calls[(calls.gnomad_freq_overall <= 0.001) | (calls.Pathogenic_Name != '.')]

# Filter for SegDups and Centromere/Telomere <= 50% overlap
# Do not apply this to Pathogenic CNVs
calls = calls[(calls["%SD"] <= 0.5) | (calls.Pathogenic_Name != '.')]
calls = calls[(calls["%CenTel"] <= 0.5) | (calls.Pathogenic_Name != '.')]

# Do some reformatting to match other annotations
calls.Type = calls.Type.apply(lambda row: row.upper())

# Add a variant ID
calls['variant_id'] = calls.Chromosome+'_'+calls.Start.astype(str)+'_'+calls.End.astype(str)+'_'+calls.Type+'_'+calls.PatientID

# Save to file
calls.to_csv(OUTPUT, sep = '\t', index = False)

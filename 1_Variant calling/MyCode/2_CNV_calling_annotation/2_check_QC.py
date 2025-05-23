import os
import pandas as pd

# Concat all QC summary files together

# Input and output files
pcnv_output="/path/to/PennCNV/output/directory" # Use the OUTPUT_DIR from script 1_run_PennCNV.sh
output_qc_file="/path/to/output/qc/file.txt"

qc_sum_path=f'{pcnv_output}/QC/qc_sum'

# Concat QC summary files across all batches
auto_qc=pd.DataFrame()
x_qc=pd.DataFrame()
for BATCH in ['OMNI', 'GSA']
	for chrom in ['autosome', 'chrX']:
		fdf=pd.read_csv(f'{qc_sum_path}/{chrom}_{BATCH}.qcsum', sep='\t')
		if chrom=='autosome':
			auto_qc=pd.concat([auto_qc, fdf])
		else:
			x_qc=pd.concat([x_qc, fdf])

# Helper functions
def auto_pass(row):
	lrr = row['LRR_SD']
	cnv = row['NumCNV']
	wf = row['WF']

	if lrr <= 0.35 and cnv <= 30 and wf <= 0.03 and wf >= -0.03:
		return True

	return False

def X_pass(row):
	lrr = row['LRR_XSD']
	cnv = row['NumCNV']

	if lrr <= 0.35 and cnv <= 30:
		return True

	return False

# Note which samples passed QC
auto_qc['Pass'] = auto_qc.apply(auto_pass, axis = 1)
x_qc['X_Pass'] = x_qc.apply(X_pass, axis = 1)

# Combine QC files
df = pd.merge(auto_qc, x_qc, how = 'outer', on = 'File')

# Save
df.to_csv(output_qc_file, sep='\t', index=False)

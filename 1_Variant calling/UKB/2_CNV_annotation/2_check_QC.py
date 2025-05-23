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
for (root,dirs,files) in os.walk(qc_sum_path, topdown=True):
	for f in files:
		fdf=pd.read_csv(f'{qc_sum_path}/{f}', sep='\t')
		if 'autosome' in f:
			auto_qc=pd.concat([auto_qc, fdf])
		elif 'chrX' in f:
			x_qc=pd.concat([x_qc, fdf])

# Note which samples passed QC
auto_qc['Pass_auto']=False
auto_qc.loc[((auto_qc.LRR_SD<=0.35) & (auto_qc.NumCNV<=30) & (auto_qc.WF<=0.03) & (auto_qc.WF>=-0.03)), 'Pass_auto']=True

x_qc['Pass_chrX']=False
x_qc.loc[((x_qc.LRR_XSD<=0.35) & (x_qc.NumCNV<=30)), 'Pass_chrX']=True

# Combine QC files
df = pd.merge(auto_qc, x_qc, how = 'outer', on = 'File')

# Save
df.to_csv(output_qc_file, sep='\t', index=False)

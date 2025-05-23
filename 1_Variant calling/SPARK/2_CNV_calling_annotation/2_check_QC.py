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
auto_qc.loc[((auto_qc.LRR_SD<=0.35) & (auto_qc.WF<=0.05) & (auto_qc.WF>=-0.05) & (auto_qc.BAF_drift<=0.01)), 'Pass_auto']=True

x_qc['Pass_chrX']=False
x_qc.loc[(x_qc.LRR_XSD<=0.35), 'Pass_chrX']=True

# Combine QC files
df = pd.merge(auto_qc, x_qc, how = 'outer', on = 'File')

# Add sample IDs and clean up
df['Sample']=df.File.str.split('/', expand=True)[8].str.split('.', expand=True)[0]
df['File']=df.File.str.split('gunzip -c ', expand=True)[1].str.split(' | ', expand=True)[0]
df=df[['Sample', 'Pass_auto', 'Pass_chrX', 'File',
	'LRR_mean', 'LRR_median', 'LRR_SD', 'BAF_mean', 'BAF_median', 'BAF_SD', 'BAF_drift', 'WF', 'NumCNV_x',
	'LRR_Xmean', 'LRR_Xmedian', 'LRR_XSD', 'BAF_Xhet', 'NumCNV_y']]
df.columns=['Sample', 'Pass_auto', 'Pass_chrX', 'File',
		'LRR_mean', 'LRR_median', 'LRR_SD', 'BAF_mean', 'BAF_median', 'BAF_SD', 'BAF_drift', 'WF', 'NumCNV_auto',
		'LRR_Xmean', 'LRR_Xmedian', 'LRR_XSD', 'BAF_Xhet', 'NumCNV_chrX']

# Save
df.to_csv(output_qc_file, sep='\t', index=False)

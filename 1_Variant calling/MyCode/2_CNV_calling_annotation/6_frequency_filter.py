import pandas as pd

# Filter calls for those in known dosage-sensitive regions or present at less than 0.1% frequency in a control cohort

# Input and output files
QC_file="/path/to/qc/file.txt" # Use the output of script 2_check_QC.py
INPUT_DIR="/path/to/input/bed/files/" # Use the output directory from script 5_frequency_pathogenicity_anno.sh
output_file="/path/to/output/file.bed"

# Number of control samples
ncontrol=18572

# Get the number of MyCode participants who passed QC
qc_pass = pd.read_csv(QC_FILE, sep = '\t')
pass_num = qc_pass[(qc_pass.Pass_auto==True) & (qc_pass.Pass_chrX==True)].shape[0]
# Save memory
del qc_pass

# Deletions
dels = pd.read_csv(f'{input_path}/intracohort_anno_dels.bed', sep = '\t',
					names = ['Chr', 'Start', 'End', 'Type', 'Zygosity', 'Length', 'NumSNP', 'Sample', 'StartSNP', 'EndSNP', 'Merge', 'Control_num', 'MyCode_num'])
dels['Control_freq'] = dels.Control_num / ncontrol
dels['MyCode_freq'] = dels.MyCode_num / pass_num
dels = dels[dels.Control_freq <= 0.001]
dels['variant_id'] = dels.Chr + '_' + dels.Start.map(str) + '_' +  dels.End.map(str) + '_' + dels.Type + '_' + dels.Sample.map(str)

# Duplications
dups = pd.read_csv(f'{input_path}/intracohort_anno_dups.bed', sep = '\t',
					names = ['Chr', 'Start', 'End', 'Type', 'Zygosity', 'Length', 'NumSNP', 'Sample', 'StartSNP', 'EndSNP', 'Merge', 'Control_num', 'MyCode_num'])
dups['Control_freq'] = dups.Control_num / ncontrol
dups['MyCode_freq'] = dups.MyCode_num / pass_num
dups = dups[dups.Control_freq <= 0.001]
dups['variant_id'] = dups.Chr + '_' + dups.Start.map(str) + '_' +  dups.End.map(str) + '_' + dups.Type + '_' + dups.Sample.map(str)

# Pathogenic CNVs
# Deletions
patho_del = pd.read_csv(f'{input_path}/pathogenic_anno_dels.bed', sep = '\t',
	names = ['Chr', 'Start', 'End', 'Type', 'Zygosity', 'Length', 'NumSNP', 'Sample', 'StartSNP', 'EndSNP', 'Merge', 'Control_num', 'MyCode_num', 'Patho_chr', 'Patho_Start', 'Patho_End', 'Patho_Name'])
patho_del['Control_freq'] = patho_del.Control_num / ncontrol
patho_del['MyCode_freq'] = patho_del.MyCode_num / pass_num
patho_del['variant_id'] = patho_del.Chr + '_' + patho_del.Start.map(str) + '_' +  patho_del.End.map(str) + '_' + patho_del.Type +  '_' + patho_del.Sample.map(str)

# Duplications
patho_dup = pd.read_csv(f'{input_path}/pathogenic_anno_dups.bed', sep = '\t',
	names = ['Chr', 'Start', 'End', 'Type', 'Zygosity', 'Length', 'NumSNP', 'Sample', 'StartSNP', 'EndSNP', 'Merge', 'Control_num', 'MyCode_num', 'Patho_chr', 'Patho_Start', 'Patho_End', 'Patho_Name'])
patho_dup['Control_freq'] = patho_dup.Control_num / ncontrol
patho_dup['MyCode_freq'] = patho_dup.MyCode_num / pass_num
patho_dup['variant_id'] = patho_dup.Chr + '_' + patho_dup.Start.map(str) + '_' +  patho_dup.End.map(str) + '_' + patho_dup.Type +  '_' + patho_dup.Sample.map(str)

# Merge all calls
del_calls = dels.merge(patho_del, how = 'outer',
	on = ['Chr', 'Start', 'End', 'Type', 'Zygosity', 'Length', 'NumSNP', 'Sample', 'StartSNP', 'EndSNP', 'Merge', 'Control_num', 'Control_freq', 'MyCode_num', 'MyCode_freq', 'variant_id'],
	suffixes = ('', "_Pathogenic"))
dup_calls = dups.merge(patho_dup, how = 'outer',
	on = ['Chr', 'Start', 'End', 'Type', 'Zygosity', 'Length', 'NumSNP', 'Sample', 'StartSNP', 'EndSNP', 'Merge', 'Control_num', 'Control_freq', 'MyCode_num', 'MyCode_freq', 'variant_id'],
	 suffixes = ('', "_Pathogenic"))
calls = pd.concat([del_calls, dup_calls])

# Remove extra columns
calls = calls[['Chr', 'Start', 'End', 'Type', 'Zygosity', 'Length', 'NumSNP', 'Sample', 'StartSNP', 'EndSNP', 'Merge',
				'Control_num', 'Control_freq', 'MyCode_num', 'MyCode_freq', 'variant_id', 'Pathogenic_Name']]

# Replace NaN with '.' for missing data
calls['Pathogenic_Name']=calls['Pathogenic_Name'].fillna('.')

# Sort
calls.sort_values(by = ['Chr', 'Start', 'End', 'Type', 'Sample'], inplace = True)

# Save to file
calls.to_csv(output_file, sep = '\t', index = False)

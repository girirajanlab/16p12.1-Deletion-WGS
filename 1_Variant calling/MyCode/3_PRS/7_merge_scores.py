import pandas as pd

# Input and output files
INPUT_PATH="/path/to/input/PRS/files" # Use the OUTPUT_PATH from script 6_run_ldpred2.sh
OUTPUT_FILE="/path/to/output/file.csv"

df=pd.DataFrame({'IID':[], 'FID':[]})
for pheno in ['autism', 'intelligence', 'educational_attainment', 'schizophrenia']:
	indf=pd.read_csv(f'{INPUT_PATH}/{pheno}_PRS.csv')
	df=pd.merge(df, indf, on=['IID', 'FID'], how='outer')

# Save
df.to_csv(OUTPUT_FILE, index=False)

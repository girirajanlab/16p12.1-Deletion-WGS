import pandas as pd

# Create a table of the burden of different variant classes

# Input and output files
input_file='/path/to/input/file.csv' # Use the output of script 7_loeuf_scores.py here
output_file='/path/to/output/table.csv' # Note that the output of this script will be the final coding SNV annotations

# Load files
df = pd.read_csv(input_file)

# Count variant types
outdf=pd.DataFrame(df.Sample.value_counts())
types=['lof', 'missense', 'splice']
for ty in types:
	out_df[ty]=out_df.index.map(df[df.Mut_type==ty]['Sample'].value_counts().to_dict())

df=df[df.LOEUF<=0.35]
out_df['All_coding_SNVs_LF']=out_df.index.map(df.Sample.value_counts().to_dict())
for ty in types:
	out_df[ty+'_LF']=out_df.index.map(df[df.Mut_type==ty]['Sample'].value_counts().to_dict())
out_df.fillna(0, inplace=True)

# Update column names
outdf.columns=['All coding SNVs', 'LOF', 'Missense', 'Splice',
				'All coding SNVs (LF)', 'LOF (LF)', 'Missense (LF)', 'Splice (LF)']
out_df['Sample']=out_df.index

# Save
out_df[['Sample', 'All coding SNVs', 'LOF', 'Missense', 'Splice',
		'All coding SNVs (LF)', 'LOF (LF)', 'Missense (LF)', 'Splice (LF)']].to_csv(output_file, index=False)

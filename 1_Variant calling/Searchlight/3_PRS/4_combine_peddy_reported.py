import pandas as pd
import numpy as np

# Use peddy-calculated ancestries and self reported ancestry to identify european samples

# Input and output files
peddy_res="/path/to/peddy/results.het_check.csv" # Use the output of script 3_peddy_ancestry.sh
reported_path='/path/to/Searchlight/background_history_files'
# Phenotypic data accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
output_european_sample_file="/path/to/output/file.fam"

# Load files
pred=pd.read_csv(peddy_res)

cols=['race_african_amer', 'race_asian', 'race_native_amer', 'race_native_hawaiian', 'race_other', 'race_white', 'race_more_than_one']
suffixes=['', '_adult']
rep=pd.DataFrame()
for suffix in suffixes:
	df=pd.read_csv(f'{reported_path}/svip_background_history{suffix}.csv')
	df=df[['individual']+[f'svip_background_history{suffix}.{i}' for i in cols[0:6]]]
	df.columns=['individual']+cols[0:6]
	rep=pd.concat([rep, df])

# Remove any samples that appear more than once
rep.drop_duplicates(subset='individual', inplace=True)
rep.fillna(False, inplace=True)

# Parse reported ancestry
rep['multiple']=rep[cols[0:6]].sum(axis=1)>1

rep['anc']='unknown'
rep.loc[rep.race_african_amer, 'anc']='AFR'
rep.loc[rep.race_asian, 'anc']='EAS'
rep.loc[rep.race_native_amer, 'anc']='NAMR'
rep.loc[rep.race_native_hawaiian, 'anc']='NAH'
rep.loc[rep.race_other, 'anc']='other'
rep.loc[rep.race_white, 'anc']='EUR'
rep.loc[rep.multiple, 'anc']='multiple'

rep.individual=rep.individual.str.replace('-', '.')

# Annotate self reported ancestry to predicted ancestry
pred['Sample']=pred.sample_id.str.split('_', expand=True)[1]
pred['rep_anc']=pred.Sample.map(dict(zip(rep.individual.to_list(), rep.anc.to_list())))

# Restrict to European samples
pred.fillna('unknown', inplace=True)
pred['comb_anc']=pred.rep_anc+'.'+pred['ancestry-prediction']
pred=pred[pred.comb_anc.isin(['EUR.EUR', 'unknown.EUR', 'EUR.UNKNOWN'])]

# Save EUR samples
pred[['fam', 'Sample']].to_csv(output_european_sample_file, sep=' ', index=False, header=False)


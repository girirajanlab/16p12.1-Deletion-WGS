import pandas as pd
import numpy as np

# Use peddy-calculated ancestries and self reported ancestry to identify european samples

# Input and output files
peddy_res="/path/to/peddy/results.het_check.csv" # Use the output of script 3_peddy_ancestry.sh
reported='/path/to/individuals_registration.csv'
output_european_sample_file="/path/to/output/file.fam"

# Load files
pred=pd.read_csv(peddy_res)
rep=pd.read_csv(reported)

# Parse reported ancestry
rep=rep[['subject_sp_id', 'family_sf_id',
	'race_asian', 'race_african_amer', 'race_native_amer', 'race_native_hawaiian', 'race_white', 'race_other', 'race_more_than_one_calc',
	'hispanic']]

rep['anc']='unknown'
rep.loc[rep.race_white==1, 'anc']='EUR'
rep.loc[rep.race_other==1, 'anc']='unknown'
rep.loc[rep.race_asian==1, 'anc']='SAS'
rep.loc[rep.race_african_amer==1, 'anc']='AFR'
rep.loc[rep.race_native_amer==1, 'anc']='NAMR'
rep.loc[rep.race_native_hawaiian==1, 'anc']='NAH'
rep.loc[rep.hispanic==1, 'anc']='hispanic'
rep.loc[rep.race_more_than_one_calc==1, 'anc']='multiple'

# Annotate self reported ancestry to predicted ancestry
pred['Sample']=pred.sample_id.str.split('_', expand=True)[1]
pred['rep_anc']=pred.Sample.map(dict(zip(rep.IID.to_list(), rep.anc.to_list())))

# Restrict to European samples
pred.fillna('unknown', inplace=True)
pred['comb_anc']=pred.rep_anc+'.'+pred['ancestry-prediction']
pred=pred[pred.comb_anc.isin(['EUR.EUR', 'unknown.EUR', 'EUR.UNKNOWN'])]

# Save EUR samples
pred[['fam', 'Sample']].to_csv(output_european_sample_file, sep=' ', index=False, header=False)


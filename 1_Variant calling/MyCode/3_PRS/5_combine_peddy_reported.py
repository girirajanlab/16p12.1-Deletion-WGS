import pandas as pd
import numpy as np

# Use peddy-calculated ancestries and self reported ancestry to identify european samples

# Input and output files
peddy_res="/path/to/peddy/results.het_check.csv" # Use the output of script 4_peddy_ancestry.sh
reported_ancestry="/path/to/reported/ancestry.csv" # Reported ancestry data was provided by MyCode
output_european_sample_file="/path/to/output/file.fam"

# Load files
pred=pd.read_csv(peddy_res)
rep=pd.read_csv(reported_ancestry)

# Parse reported ancestry
rep=rep[['SEQN_ID', 'PT_RACE_1', 'PT_ETHNICITY']]

rep['anc']='unknown'
rep.loc[rep.PT_RACE_1=='White', 'anc']='EUR'
rep.loc[rep.PT_RACE_1=='Black Or African American', 'anc']='AFR'
rep.loc[rep.PT_ETHNICITY=='Hispanic or Latino', 'anc']='hispanic'

# Annotate self reported ancestry to predicted ancestry
pred['Sample']=pred.sample_id.str.split('_GHS', expand=True)[1]
pred['rep_anc']=pred.Sample.map(dict(zip(rep.SEQN_ID.to_list(), rep.anc.to_list())))

# Restrict to European samples
pred.fillna('unknown', inplace=True)
pred['comb_anc']=pred.rep_anc+'.'+pred['ancestry-prediction']
pred=pred[pred.comb_anc.isin(['EUR.EUR', 'unknown.EUR', 'EUR.UNKNOWN'])]

# Save EUR samples
pred[['fam', 'Sample']].to_csv(output_european_sample_file, sep=' ', index=False, header=False)


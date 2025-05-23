import pandas as pd
import numpy as np

# Use peddy-calculated ancestries and self reported ancestry to identify european samples

# Input and output files
peddy_res="/path/to/peddy/results.het_check.csv" # Use the output of script 3_peddy_ancestry.sh
reported='/path/to/UKB/self/reported/race/data.csv' # A file that contains data fields 'eid', '21000-0.0', '21000-1.0', and '21000-2.0'
output_european_sample_file="/path/to/output/file.fam"

# Load files
pred=pd.read_csv(peddy_res)
rep=pd.read_csv(reported)

# Parse reported ancestry
rep['combined']=rep[['21000-0.0', '21000-1.0', '21000-2.0']].apply(lambda row: ';'.join(list(set([str(int(i)) for i in row.to_list() if i==i]))), axis=1)
rep['first_digit']=rep.combined.str[0]

rep['anc']='unknown'
rep.loc[rep.first_digit=='1', 'anc']='EUR'
rep.loc[rep.first_digit=='2', 'anc']='multiple'
rep.loc[rep.first_digit=='3', 'anc']='SAS'
rep.loc[rep.first_digit=='4', 'anc']='AFR'
rep.loc[rep.first_digit=='5', 'anc']='EAS'

# Annotate self reported ancestry to predicted ancestry
pred['rep_anc']=pred.Sample.map(dict(zip(rep.eid.to_list(), rep.anc.to_list())))

# Restrict to European samples
pred.fillna('unknown', inplace=True)
pred['comb_anc']=pred.rep_anc+'.'+pred['ancestry-prediction']
pred=pred[pred.comb_anc.isin(['EUR.EUR', 'unknown.EUR', 'EUR.UNKNOWN'])]

# Save EUR samples
pred['fam']=pred.sample_id.str.split('_', expand=True)[0]
pred[['fam', 'Sample']].to_csv(output_european_sample_file, sep=' ', index=False, header=False)


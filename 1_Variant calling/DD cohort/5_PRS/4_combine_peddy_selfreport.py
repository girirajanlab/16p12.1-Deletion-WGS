import pandas as pd
import numpy as np

# Use peddy-calculated ancestries and self reported ancestry to identify european samples

# Input and output files
peddy_res="/path/to/peddy/results.het_check.csv" # Use the output of script 3_peddy_ancestry.sh
selfreport='/path/to/self_reported_ancestry.csv'
output_european_sample_file="/path/to/output/file.fam"

# Load files
pred=pd.read_csv(peddy_res)
rep=pd.read_csv(selfreport)

# Parse self reported ancestry
mapping={np.nan:'unknown', '.':'unknown', 'Caucasian':'EUR', 'British':'EUR', 'White':'EUR', 'Australian':'EUR', 'European':'EUR',
	'Indian':'SAS', 'Aboriginal Australian':'AAUS', 'White/Jewish':'EUR', 'African':'AFR', 'English':'EUR', 'Chinese/European':'EAS/EUR',
	'Chinese':'EAS', 'Australian/Italian':'EUR', 'Asian':'EAS', 'Caucasian/American Indian':'EUR/NAMR', 'Caucasian/Lebanse':'WANA', 'Lebanese':'WANA'}
rep['sr_anc']=rep.self_reported_ancestry.map(mapping)

# Annotate self reported ancestry to predicted ancestry
pred['Sample']=pred.sample_id.str.split('_', expand=True)[1]
pred['sr_anc']=pred.Sample.map(dict(zip(rep.IID.to_list(), rep.sr_anc.to_list())))

# Restrict to European samples
pred.fillna('unknown', inplace=True)
pred['comb_anc']=pred.sr_anc+'.'+pred['ancestry-prediction']
pred=pred[pred.comb_anc.isin(['EUR.EUR', 'unknown.EUR', 'EUR.UNKNOWN'])]

# Save EUR samples
pred[['fam', 'Sample']].to_csv(output_european_sample_file, sep=' ', index=False, header=False)


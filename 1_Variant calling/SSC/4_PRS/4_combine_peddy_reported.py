import pandas as pd
import numpy as np

# Use peddy-calculated ancestries and self reported ancestry to identify european samples

# Input and output files
peddy_res="/path/to/peddy/results.het_check.csv" # Use the output of script 3_peddy_ancestry.sh
reported_path='/path/to/SSC/core/descriptive/files'
# Phenotypic data accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
output_european_sample_file="/path/to/output/file.fam"

# Load files
pred=pd.read_csv(peddy_res)

types=['proband', 'mother', 'father', 'sibling', 'sibling2', 'twin']
rep=pd.DataFrame()
for t in types:
	df=pd.read_csv(f'{reported_path}/{t}_ssc_core_descriptive.csv')
	df=df[['individual', 'race', 'ethnicity']]
	rep=pd.concat([rep, df])

# Parse reported ancestry
mapping={'white':'EUR', 'more-than-one-race':'multiple', 'other':'other', 'asian':'SAS', 'african-amer':'AFR', 'not-specified':'unknown', 'native-american':'NAMR', 'native-hawaiian':'NAH', np.nan:'unknown'}
rep['anc']=rep.race.map(mapping)
rep.loc[rep.ethnicity=='hispanic', 'anc']='hispanic'

# Annotate self reported ancestry to predicted ancestry
pred['Sample']=pred.sample_id.str.split('_', expand=True)[1]
pred['rep_anc']=pred.Sample.map(dict(zip(rep.individual.to_list(), rep.anc.to_list())))

# Restrict to European samples
pred.fillna('unknown', inplace=True)
pred['comb_anc']=pred.rep_anc+'.'+pred['ancestry-prediction']
pred=pred[pred.comb_anc.isin(['EUR.EUR', 'unknown.EUR', 'EUR.UNKNOWN'])]

# Save EUR samples
pred[['fam', 'Sample']].to_csv(output_european_sample_file, sep=' ', index=False, header=False)


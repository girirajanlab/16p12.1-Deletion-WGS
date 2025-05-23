import pandas as pd

# Gather SRS data from Searchlight probands

# Input and output files
SRS_PARENT="/path/to/Simons_Searchlight_Phase1_16p11.2_Dataset_v11.0/srs_parent.csv"
SUBJECTS="/path/to/Simons_Searchlight_Phase1_16p11.2_Dataset_v11.0/svip_subjects.csv"

OUTPUT_FILE="/path/to/output/file/for/SRS/scores.csv"

# Load SRS data
searchlight=pd.read_csv(SRS_PARENT)
searchlight=searchlight[['individual', 'family_type', 'srs_parent.total']]
searchlight.columns=['Sample', 'family_type', 'SRS Raw Score']
searchlight=searchlight[searchlight.family_type.isin(['16p-deletion', '16p-duplication'])]
searchlight=searchlight[~searchlight['SRS Raw Score'].isnull()]

# Restrict to probands
searchlight_samples=pd.read_csv(SUBJECTS)
searchlight_samples=searchlight_samples[searchlight_samples.relationship_to_iip=='Initially identified proband']

searchlight=searchlight[searchlight.Sample.isin(searchlight_samples.sfari_id.to_list())]

# Save
searchlight.to_csv(OUTPUT_FILE, index=False)
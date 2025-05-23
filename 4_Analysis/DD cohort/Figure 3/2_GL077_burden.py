import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42

# Plot the burden in carriers of 16p12.1 deletion in family GL077

# Input and Output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
OUTPUT_FIG="/path/to/output/boxplot/figure.pdf" # This plot will be the lineplot in Fig 3B left

# Load data
df=pd.read_csv(TABS1A)
df=df[df.Sample.isin(['GMC_077', 'MC_077', 'P1C_077'])]
df=df[['Sample', 'All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs']]

# Convert to long
long_df=df.melt(id_vars='Sample')

sns.lineplot(data=long_df, x='Sample', y='value', hue='variable', palette=['#E50E2E', '#93B8A3', '#4B4D7F', '#33977C'], hue_order=['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs'])
sns.scatterplot(data=long_df, x='Sample', y='value', hue='variable', s=100, palette=['#E50E2E', '#93B8A3', '#4B4D7F', '#33977C'], hue_order=['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs'], marker='o', facecolor='white', legend=False)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig(OUTPUT_FIG)
import pandas as pd
import numpy as np
import statsmodels.stats.power as smp

import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['pdf.fonttype']=42

# Calculate the power for differences in burden between probands and parents

# Input and output files
SUMM_STAT="Power_summary_statistics.xlsx" # Summary statistics gathered from published data and in-house analysis
OUTPUT_STATS="/path/to/output/power/statistics.csv"
OUTPUT_FIG="/path/to/output/figure.pdf" # The plot presented in Fig S2B

# Calculate power
summ_stat=pd.read_excel(SUMM_STAT)
summ_stat=summ_stat[~summ_stat.Source.isnull()]

# Define the Cohen's D effect size
def cohens_d(g1m, g1s, g1n, g2m, g2s, g2n):
    dof = g1n + g2n - 2
    return (g1m - g2m) / np.sqrt(((g1n-1)*g1s ** 2 + (g2n-1)*g2s ** 2) / dof)

summ_stat["Cohen's D"]=summ_stat.apply(lambda x: cohens_d(x['Group 1 Mean'], x['Group 1 SD'], x['Group 1 N'], x['Group 2 Mean'], x['Group 2 SD'], x['Group 2 N']), axis=1)
print(summ_stat)

vars=summ_stat.Variant.to_list()
power_out=[]
for v in vars:
	d=summ_stat[summ_stat.Variant==v]["Cohen's D"].to_list()[0]
	alpha=0.05
	alternative='larger'
	for n in range(10, 160, 10):
		power=smp.TTestPower().power(effect_size=d, nobs=n, alpha=alpha, alternative=alternative)
		power_out.append([v, n, d, alpha, alternative, power])
power_df=pd.DataFrame(power_out, columns=['Variant', 'N', "Cohen's D", "Alpha", "Alternative", "Power"])

# Save
power_df.to_csv(OUTPUT_STATS, index=False)

# Plot power

# Add some dummy variables for plot formatting
power_df['style']=1

sns.lineplot(data=power_df, x="N", y="Power", hue="Variant")
sns.scatterplot(data=power_df, x="N", y="Power", hue="Variant", legend=False)

# Add a line at 54
plt.plot([54, 54], [0, 2], color='k', ls='--', zorder=0)
plt.ylim(0, 1.1)

# Add another line at 80% power
plt.plot([0, 160], [0.8, 0.8], color='k', ls='--', zorder=0)
plt.xlim(0, 160)

plt.xlabel('Sample size')

plt.savefig(OUTPUT_FIG)
This README describes the scripts used for analysis of samples in the DD cohort

# Scripts
This lists all of the scripts in the directory and subdirectories and a description
- Figure 2/
	- _1_bar_graphs.py_: Create bar graphs showcasing phenotypes in the DD cohort (Fig 2A-B)
	- _2_proband_v_parent_shift.py_: Create KDE plots showing the differences in quatitative phenotypes between DD cohort probands and their parents (Fig 2C top)
	- _3_DD_vs_SSC_Searchlight.py_: Create KDE plots showing the differences in quantiative phenotypes between DD cohort probands and autism probands from SSC and Searchlight (Fig 2C middle)
	- _4_bmi_hc.py_: Create KDE plots showing BMI and head circumference distributions in DD cohort probands (Fig 2C bottom)
	- _5_milestone_boxplots.py_: Create boxplots showing differences in milestone achievement from carrier and non-carrier children from the DD cohort (Fig 2D)
- Figure 3/
	- _1_paired_ttest.py_: Compare proband and carrier and noncarrier parent burden using paired t-tests (Fig 3A and S2C)
	- _2_GL_077_burden.py_: Examine the burden in GL_077 (Fig 3B)
	- _3_brian_celltype_enrichment.py_: Calculate enrichment of secondary variants in genes preferentially expressed in different neuronal subtypes
	- _4_WGCNA_enrichment.py_: Calculate enrichment of secondary variant genes in genes coexpressed with 16p12.1 deletion genes in NPCs derived from DD cohort families
- Figure 4/
	- _1_16p12_connectors.py_: Identify genes along shortest paths between 16p12.1 genes and DD cohort proband secondary variant genes and identify GO enrichments on unique connector genes with gProfiler
	- _2_find_hubs.py_: Find network hub genes per DD cohort proband and create GML files for plotting
- Figure 5/
	- _1_regressions.py_: Perform linear and logistic regressions assessing effects of secondary variants on phenotpic domains and quantitative phenotypes in DD cohort probands
	- _2_gProfiler.py_: Perform GO and KEGG enrichments on secondary variant genes from probands with each phenotypic domain
	- _3_geneset_ttests.py_: Assess enrichment of secondary variant genes from probands with each phenotypic domain from curated disease gene sets
- Figure 6/
	- _1_adult_assoctaions.py_: Perform t-tests comparing variant burden in DD cohort 16p12.1 deletion carrier adults with and without specific phenotypes
	- _2_child_associations.py_: Perform t-tests comparing variant burden in DD cohort 16p12.1 deletion carrier children with and without specific phenotypes
	- _3_phenotype_prevalence.py_: Calulate the prevalence of specific phenotypes in DD cohort adult and child 16p12.1 deletion carriers

- Figure S2/
	- _1_identify_pathogenic_SNVs.py_: Identify pathogenic SNVs in DD cohort probands
	- _2_updated_pathogenic_SNVs.py_: Compile pathogenic variants in DD cohort probands
- Figure S3/
	- _1_geneset_enrichment.py_: Calculate enrichment of DD cohort proband secondary variants for gene from curated disease gene sets
	- _2_brainspan_enrichment.py_: Calculate spatial and temporal enrichment of secondary variant genes in the developing human brain
- Figure S4/
	- _1_burden_ttest.py_: Compare the burden of each variant class in probands with and without each phenotypic domain
	- _2_quantitative_phenotype_correlation.py_: Calculate Pearson correlations between variant burden and quantitative phenotypes
- Figure S5/
	- _1_DD_vs_Estonia.py_: Compare burden in 16p12.1 deletion probands and parents from the DD cohort to carriers in the Estonian Biobank
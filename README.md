# 16p12.1-Deletion-WGS
Scripts for:
Jensen, M., Smolen, C., Tyryshkina, A., Pizzo, L., Banerjee, D., Oetjens, M., Shimelis, H., Taylor, C. M., Pounraja, V. K., Song, H., Rohan, L., Huber, E., El Khattabi, L., van de Laar, I., Tadros, R., Bezzina, C., van Slegtenhorst, M., Kammeraad, J., Prontera, P., … Girirajan, S. (2024). Genetic modifiers and ascertainment drive variable expressivity of complex disorders. _medRxiv: The Preprint Server for Health Sciences_, 2024.08.27.24312158. https://doi.org/10.1101/2024.08.27.24312158

# Repository Organization
This repository contains scripts for calling genetic variants and performing analyses in the above manuscript. Files are organized into four main directories, each with additional READMEs. There are too many scripts to detail here, so all scripts will be detailed in the READMEs within each subdirectory. Some examples will be shown here. The four main directories are:
  ## 1_Variant calling/
  Contains scripts for calling variants in the DD, MyCode, Simons Searchlight, Simons SSC, Simons SPARK, and UK Biobank cohorts. While variant calling is very similar across cohorts, files are separated by cohort and then by variant type, each with a separate README file. As an example, the DD cohort directory contains five sub-folders, each with its own README.
  
 - DD cohort/
   - 1_GATK/
   - 2_SNV_annotation/
   - 3_CNV_calling_annotation/
   - 4_STR_calling_annotation/
   - 5_PRS/

  ## 2_Analysis preparation/
  Contains scripts related to processing of other data relevant for analyses in this manuscript.
  - Brain_Celltype_Expression/
    - _1_parse_Allen_data.py_: Identifying preferentially expressed genes in specific brain cell types using data from the Allen Brain Institute
  - Gene_Annotations/
    - _0_download_annotations.sh_: Download gene annotations for identifying disease-associated genes
    - _1_werling_annotations.py_: Parse disease gene annotations from [Werling et al Nature Genetics 2018](https://pubmed.ncbi.nlm.nih.gov/29700473/)
    - _2_dbd_ddg2p_sfari.py_: Parse disease gene annotations from the DBD Genes Database, Gene2Phenotype (G2P) database, and SFARI Gene
    - _3_szdb_anno.py_: Parse schizophrenia annotations from SZDB using GWAS, CNVs, or Exome sequencing results
    - _4_epilepsy_anno.py_: Parse epilepsy-associated gene annotations from [Wang et al. Seizure 2017](https://pubmed.ncbi.nlm.nih.gov/28007376/)
    - _5_add_loeuf.py_: Add gnomAD LOEUF annotations to genes
  - NPC_RNAseq/
    - _1_kallisto.sh_: Trim reads and quantify transcript abundance with kallisto
    - _2_WGCNA.R_: Perform WGCNA on RNA-seq data
  - STRING_Network/
    -  _1_parse_string.py_: Annotate genes and updated interaction scores in STRING DB
    -  _2_string_connectivity.py_: Quantify connectivity of genes in annotated STRING network

  ## 3_Data preparation/
  Describes the integration of genetic and phenotypic data for the DD, MyCode, Simons Searchlight, Simons SSC, Simons SPARK, and UK Biobank cohorts. Scripts are separated by cohort. As an example, the DD cohort directory is organized as:
  - DD cohort/
    - _1_make_table_s1a.py_: Integration of genetic and phenotypic data for samples in the DD cohort. The output of this script is available as Table S1A.
    - _2_make_genelists.py_: Generation of lists of genes with secondary variants for DD cohort probands, both collectively and in groups of probands with specific phenotypic domains.
    - _3_genelist_by_proband.py_: Generations of lists of secondary variant genes for each DD cohort proband.

  ## 4_Analysis/
  Contains scripts for analyses described in the manuscript. Scripts are separated by cohort and figure. There are separate folders for the DD, MyCode, Simons Searchlight, Simons SSC, Simons SPARK, and UK Biobank cohorts and additional folders for analyses that span multiple cohorts (MultiCohort) and a separate folder for a power analysis (Power analysis). There is a separate README for each cohort, describing all of the scripts and analyses for each cohort. As an example, the Simons Searchlight directory is organized as:
  - Searchlight/
    - Figure 4/
      - _1_16p11.2_connectors.py_: Identify shortest paths between 16p11.2 region genes and secondary variant genes in 16p11.2 deletion probands.
    - Figure 7/
      - _1_regressions.py_: Perform linear regressions assessing effects of secondary variants on quantitative phenotypes in 16p11.2 CNV probands
      - _2_gProfiler.py_: Perform GO and KEGG enrichments on secondary variant genes from 16p11.2 CNV probands
    - Figure S6/
      - _1_variant_correlations.py_: Calculate Pearson correlations for variant burden and quantitative phenotypes in 16p11.2 CNV probands

**Note**:  We note that data from All of Us was analyzed in the [All of Us Researcher Workbench](https://researchallofus.org/data-tools/workbench/). Scripts for All of Us are very similar to those for UK Biobank. We intend to make our AoU workspace a community workspace soon to make the scripts accessible to all registered AoU researchers.

# Acknowledgements
We thank the participants and investigators in the the DD, MyCode, Simons Searchlight, Simons SSC, Simons SPARK, UK Biobank, and All of US research studies who made this work possible. This research has been conducted using the UK Biobank Resource under Application Number 45023 and All of Us Controlled Tier Dataset v8. We also thank the National Institutes of Health’s All of Us Research Program for making available the participant data examined in this study.

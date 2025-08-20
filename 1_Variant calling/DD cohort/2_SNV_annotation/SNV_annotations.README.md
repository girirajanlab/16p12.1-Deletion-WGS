This README describes the scripts used for SNV annotations in the DD cohort

# Scripts
This lists all of the scripts in the directory and subdirectories and a description
- _1_filter_quality.sh_: Perform QC on variants from GATK
- _2_annovar.sh_: Annotate variants with functional consequences using ANNOVAR
- coding_annotations/
	- _3_filter_exonic.sh_: Filter variants for those with exonic and splicing functions
	- _4_annotate_gnomad.sh_: Annotate variants for gnomAD allele frequencies using VCFanno
	- _5_somalier_ancestry.sh_: Calculate ancestry for samples using somalier
	- _6_combine_somalier_selfreport.py_: Merge somalier-calculated ancestries with self-reported data
	- _7_filter_pop_frequency.sh_: Filter variants for those with ancestry-specific allele frequencies < 0.001
	- _8_annotate_cadd.sh_: Annotate variants with CADD scores using VCFanno
	- _9_annotate_ClinVar.sh_: Annotate variants with ClinVAR variants using VCFanno
	- _10_vcf2table.sh_: Merge all VCFs across samples into a single table
	- _11_filter_variant_types.py_: Filter variants for missense or splice variants with CADD>=25 or LOF variants
	- _12_intracohort_filter.py_: Apply an intracohort filter such that each variant occurs in <=10 individuals in the cohort
	- _13_annotate_gencode_genes.py_: Annotate variants with GENCODE genes
	- _14_loeuf_scores.py_: Annotate variants with LOEUF scores
	- _15_inheritance_annotations.py_: Annotate variants with inheritance
- noncoding_annotations/
	- enhancer/
		- _3_annotate_fetal_brain_enhancer.sh_: Annotate variants with whether they occur in putative fetal brain active enhancer regions
		- _4_filter_enhancers.sh_: Filter variants for those located in putative fetal brain active enhancer regions
		- _5_annotate_gnomad.sh_: Annotate variants for gnomAD allele frequencies using VCFanno
		- _6_filter_gnomad.sh_: Filter variants for those with ancestry-specific allele frequencies < 0.001
		- _7_vcf2table.sh_: Merge all VCFs across samples into a single table
		- _8_intracohort_filter.py_: Apply an intracohort filter such that each variant occurs in <=10 individuals in the cohort
	- promoter_UTR5/
		- _2.1_annovar_upstream.sh_: Annotate variants with functional consequences using ANNOVAR, with a focus on up/downstream variants
		- _3_filter_promoter.sh_: Filter variants for those within 5' UTRs or upstream of genes
		- _4_annotate_gnomad.sh_: Annotate variants for gnomAD allele frequencies using VCFanno
		- _5_filter_gnomad.sh_: Filter variants for those with ancestry-specific allele frequencies < 0.001
		- _6_vcf2table.sh_: Merge all VCFs across samples into a single table
		- _7_format_table.sh_: Clean up variant table before further processing
		- _8_intracohort_filter.py_: Apply an intracohort filter such that each variant occurs in <=10 individuals in the cohort

# Files
This lists all of the non-script annotation files in the directory and subdirectories and a description
- coding_annotations/
	- toml_files/
		- _gnomad_annotation_file.toml_: A TOML file describing annotation of gnomAD allele frequencies for VCFanno
		- _cadd_annotation_file.toml_: A TOML file describing annotation of CADD scores for VCFanno
		- _clinvar_annotation_file.toml_: A TOML file describing annotation of ClinVar alleles for VCFanno
- noncoding_annotations/
	- enhancer/
		- annotation_files/
			- _fetal_brain_enhancers_annotated.bed.gz[.tbi]_: A merge of Roadmap Epigomics HMM chromatin states 6, 7, and 12--which correspond with enhancer sites--in the fetal brain
			- toml_files/
				- _fetal_brain_enhancer.toml_: A TOML file describing annotation of putative fetal brain active enhancer regions for VCFanno
				- _gnomad_annotation_file.toml_: A TOML file describing annotation of gnomAD allele frequencies for VCFanno
	- promoter_UTR5/
		- toml_files/
			- _gnomad_annotation_file.toml_: A TOML file describing annotation of gnomAD allele frequencies for VCFanno

# Pipeline to process GWAS summary statistics
IDS, = glob_wildcards("data/raw/{id}.gwas.imputed_v3.both_sexes.tsv.bgz")
MAF = ["maf_0.01"]
PVAL = ["pval_10e-2000"]

def get_params(x):
	out = x.split("_")[1]
	return out

rule all:
	input:
		expand("data/iaf/{maf}/{pval}/{id}.gwas.imputed_v3.both_sexes.tsv", id=IDS, maf=MAF, pval=PVAL)


#This rule unzips raw bgz files
rule unzip_bgz:
	input:
		"data/raw/{id}.gwas.imputed_v3.both_sexes.tsv.bgz"
	output:
		"data/unzipped/{id}.gwas.imputed_v3.both_sexes.tsv"
	shell:
		"gunzip -c {input} > {output}"

#This rule runs the get_iaf file
rule iaf:
	input:
		ss = "data/unzipped/{id}.gwas.imputed_v3.both_sexes.tsv",
		blocks = "data/fourier_ls-all_parsed.bed",
		codes = "data/ukbb_phenotype_names.tsv"
	params:
		maf_t = lambda wildcards: get_params(wildcards.maf),
		pval_t = lambda wildcards: get_params(wildcards.pval)
	output:
		"data/iaf/{maf}/{pval}/{id}.gwas.imputed_v3.both_sexes.tsv"
	shell:
		"Rscript code/get_iaf.R {input.ss} {input.blocks} {input.codes} {params.maf_t} {params.pval_t} {output}"



	

	
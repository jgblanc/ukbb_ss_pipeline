# Pipeline to process GWAS summary statistics
IDS, = glob_wildcards("output/file_names/all/{id}.gwas.imputed_v3.both_sexes.tsv.bgz")
MAF = ["maf_0.01"]
PVAL = ["pval_1e-5", "pval_5e-8"]
THRESHOLDS = [4,5,6,7,8]

def get_params(x):
	out = x.split("_")[1]
	return out

def get_p_list(x):
  A = [str(i) for i in x]
  print(A)
  out = "-".join(A)
  print(out)
  return out


rule all:
	input:
		expand("data/snp_number/{maf}/{id}.numSNP", id=IDS, maf=MAF, pval=PVAL)


#This rule unzips raw bgz files
#rule unzip_bgz:
#	input:
#		"data/raw/{id}.gwas.imputed_v3.both_sexes.tsv.bgz"
#	output:
#		"data/unzipped/{id}.gwas.imputed_v3.both_sexes.tsv"
#	shell:
#		"gunzip -c {input} > {output}"

#This rule runs the get_iaf file
#rule iaf:
#	input:
#		ss = "data/unzipped/{id}.gwas.imputed_v3.both_sexes.tsv",
#		blocks = "data/fourier_ls-all_parsed.bed",
#		codes = "data/ukbb_phenotype_names.tsv"
#	params:
#		maf_t = lambda wildcards: get_params(wildcards.maf),
#		pval_t = lambda wildcards: get_params(wildcards.pval)
#	output:
#		"data/iaf/{maf}/{pval}/{id}.gwas.imputed_v3.both_sexes.tsv"
#	shell:
#		"Rscript code/get_iaf.R {wildcards.id} {input.blocks} {input.codes} {params.maf_t} {params.pval_t} {output}"
		
rule iaf_weighted:
	input:
#		ss = "data/unzipped/{id}.gwas.imputed_v3.both_sexes.tsv",
		blocks = "data/fourier_ls-all_parsed.bed",
		codes = "data/ukbb_phenotype_names.tsv"
	params:
		maf_t = lambda wildcards: get_params(wildcards.maf),
		pval_t = lambda wildcards: get_params(wildcards.pval)
	output:
		"data/iaf/weighted/{maf}/{pval}/{id}.gwas.imputed_v3.both_sexes.tsv"
	shell:
		"Rscript code/get_iaf_weighted.R {wildcards.id} {input.blocks} {input.codes} {params.maf_t} {params.pval_t} {output}"

rule snp_number:
	input:
#		ss = "data/unzipped/{id}.gwas.imputed_v3.both_sexes.tsv",
		blocks = "data/fourier_ls-all_parsed.bed",
		codes = "data/ukbb_phenotype_names.tsv"
	params:
		maf_t = lambda wildcards: get_params(wildcards.maf),
		pval_t = get_p_list(THRESHOLDS)
	output:
		"data/snp_number/{maf}/{id}.numSNP"
	shell:
		"Rscript code/get_snp_number.R {wildcards.id} {input.blocks} {input.codes} {params.maf_t} {params.pval_t} {output}"

	

	
# Pipeline to process GWAS summary statistics
IDS, = glob_wildcards("data/raw/{id}.gwas.imputed_v3.both_sexes.tsv.bgz")


rule all:
	input:
		expand("data/iaf/{id}.gwas.imputed_v3.both_sexes.tsv", id=IDS)


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
		blocks = "data/fourier_ls-all_parsed.bed"
	output:
		"data/iaf/{id}.gwas.imputed_v3.both_sexes.tsv"
	shell:
		"Rscript code/get_iaf.R {input.ss} {input.blocks} {output}"



	

	
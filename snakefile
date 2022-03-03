# Pipeline to process GWAS summary statistics
IDS, = glob_wildcards("data/raw/{id}.gwas.imputed_v3.both_sexes.tsv.bgz")


rule all:
	input:
		expand("data/unzipped/{id}.gwas.imputed_v3.both_sexes.tsv", id=IDS)


#This rule unzips raw bgz files
rule unzip_bgz:
	input:
		"data/raw/{id}.gwas.imputed_v3.both_sexes.tsv.bgz"
	output:
		"data/unzipped/{id}.gwas.imputed_v3.both_sexes.tsv"
	shell:
		"gunzip -c {input} > {output}"



	

	
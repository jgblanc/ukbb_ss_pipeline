
df <- fread("~/ukbb_ss_pipeline/data/ukbb_phenotype_names.tsv", header = F)
small <- sample_n(df, 20)

for (i in 1:nrow(small)) {
  
  print(i)
  name <- small[i,1]
  path <- "~/ukbb_ss_pipeline/output/file_names/small/"
  fp <- paste0(path, name, ".gwas.imputed_v3.both_sexes.tsv.bgz")
  cmd <- paste("touch", fp)
  
  system(cmd)
}

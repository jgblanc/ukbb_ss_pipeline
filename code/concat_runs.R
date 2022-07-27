library(dplyr)
library(data.table)


maf <- c("maf_0.01")
pval <- c("pval_5e-8")
names_raw <- list.files("../output/file_names/all/")
names <- rep(NA, length(names_raw))
for (i in 1:length(names_raw)) {
  names[i] <- tools::file_path_sans_ext(names_raw[i])
}
dat <- expand.grid(names, maf, pval)
colnames(dat) <- c("filename", "maf", "pval")

agg_all_data <- function(dir_path, maf, pval, filename) {
  
  df <- fread(paste0(dir_path,"/", maf, "/", pval, "/", filename))
  
  return(df)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../data/iaf/weighted/' )

fwrite(df, "../output/results/all_weighted_5e-8_0.01.txt", row.names=F,quote=F,sep="\t", col.names = T)
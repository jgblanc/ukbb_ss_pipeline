library(dplyr)
library(data.table)


maf <- c("maf_0.01")
pval <- c("pval_1e-5")
names_raw <- list.files("../data/mle/maf_0.01/pval_1e-5/")
names <- rep(NA, length(names_raw))
for (i in 1:length(names_raw)) {
  names[i] <- tools::file_path_sans_ext(names_raw[i])
}
print(head(names))
dat <- expand.grid(names, maf, pval)
colnames(dat) <- c("filename", "maf", "pval")

agg_all_data <- function(dir_path, maf, pval, filename) {
  
  df <- fread(paste0(dir_path,"/", maf, "/", pval, "/", filename, ".mle"))
  
  return(df)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../data/mle' )

fwrite(df, "../output/results/mle_lenient.txt", row.names=F,quote=F,sep="\t", col.names = T)
# This script will assign each variant to an ld block and filter to only one 
# significant SNP per block.
setwd("~/ukbb_ss_pipeline/code/")
library(data.table)
library(tidyverse)




# Read in summary statistics data
df <- fread("../data/unzipped/50_raw.gwas.imputed_v3.both_sexes.tsv")

# Remove low confidence variants
df_filtered <- df %>% filter(low_confidence_variant == FALSE)

# Remove non-significant variants
df_sig <- df_filtered %>% filter(pval < 5e-8)

# Separated column 1 into chromosome, base pair, reference allele, and alt. allele
df_split <- df_sig %>% separate(variant, c("chr", "bp", "ref", "alt"), sep = ":")

# Filtered out low minor allele frequency variants
df_af_filter <- df_split %>% filter(minor_AF > 0.01)

# Read in ld block information
ld_block <- fread("../data/fourier_ls-all_parsed.bed")

#Assigns each SNP to a specific block
assign_SNP_to_block <- function(snp_info, block) {
  
  # Get chromosome and base pair info
  CHR <- as.numeric(snp_info[,"chr"])
  BP <- as.numeric(snp_info[, "bp"])
  
  # Filter block based on SNP
  block_chr <- block %>% filter(chr == CHR)
  first_start <- as.numeric(block_chr[1,"start"])
  block_bp <- block_chr %>% filter((start < BP & stop >= BP) | BP == first_start)
  
  # Adds a column to SNP info with block number
  snp_info$block <- as.numeric(block_bp[,"block_number"])
  return(snp_info)
  
}


for (i in 1:nrow(df_af_filter)) {
  
  df_af_filter$block <- rep(0, nrow(df_af_filter))
  df_af_filter[i,] <- assign_SNP_to_block(snp_info = df_af_filter[i,], block = ld_block)
  
}

assign_SNP_to_block(snp_info = df_af_filter[2,], block = ld_block)

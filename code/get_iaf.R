# This script will assign each variant to an ld block and filter to only one 
# significant SNP per block.


args=commandArgs(TRUE)

if(length(args) != 3) {stop("Please supply 3 inputs")}

library(data.table)
library(tidyverse)
library(boot)

ss_file=args[1]
blocks_file=args[2]
output_file=args[3]

# Read in summary statistics data
df <- fread(ss_file)

# Remove low confidence variants
df_filtered <- df %>% filter(low_confidence_variant == FALSE)

# Remove non-significant variants
df_sig <- df_filtered %>% filter(pval < 5e-8)

# Separated column 1 into chromosome, base pair, reference allele, and alt. allele
df_split <- df_sig %>% separate(variant, c("chr", "bp", "ref", "alt"), sep = ":")

# Filtered out low minor allele frequency variants
df_af_filter <- df_split %>% filter(minor_AF > 0.01)

# Read in ld block information
ld_block <- fread(blocks_file)

assign_SNP_to_block <- function(CHR, BP, block = ld_block) {
  
  # Filter blocks based on snp
  block_chr <- block %>% filter(chr == CHR)
  first_start <- as.numeric(block_chr[1, "start"])
  block_bp <- block_chr %>% filter( (start < BP & stop >= BP) | BP == first_start) 
  
  # Assign
  block_num <- as.numeric(block_bp[,"block_number"])
  return(block_num)
}


df_af_filter <- df_af_filter[1:1000,] #Delete this line for final implementation
# Add block info - takes a while 
df_blocks <- df_af_filter %>% 
  mutate(block = apply(., MARGIN = 1, FUN = function(params)assign_SNP_to_block(as.numeric(params[1]), as.numeric(params[2]))))

# Pick the minimum p-value per block
df_minP <- df_blocks %>%
  group_by(block) %>%
  slice_min(pval, with_ties = F)

# Find the alternate allele frequency
df_alt_freq <- df_minP %>%
  mutate(alt_AF = case_when(alt == minor_allele ~ minor_AF, alt != minor_allele ~ (1 - minor_AF)))

# Find the increasing allele frequency
df_inc_freq <- df_alt_freq %>%
  mutate(inc_AF = case_when((beta >= 0) ~ alt_AF, (beta < 0) ~ (1 - alt_AF)))

# Find average of the increasing allele frequency
avg_inc_AF <- mean(df_inc_freq$inc_AF)

# Boot strap resampling 
boot_mean <- function(original_vector, resample_vector) {
  mean(original_vector[resample_vector])
}

# R is number of replications
mean_results <- boot(df_inc_freq$inc_AF, boot_mean, R = 5000)
ci <- boot.ci(boot.out = mean_results, type = "basic")

# Create the confidence interval
lb <- ci$basic[4]
ub <- ci$basic[5]

# Create output data for trait
df_out <- as.data.frame(matrix(c(avg_inc_AF, lb, ub), nrow = 1))
colnames(df_out) <- c("mean_iaf", "lower_ci", "upper_ci")


# Save output
fwrite(df_out, output_file, row.names = F, col.names = T, sep = "\t")

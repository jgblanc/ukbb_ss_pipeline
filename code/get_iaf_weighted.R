# This script will assign each variant to an ld block and filter to only one 
# significant SNP per block.


args=commandArgs(TRUE)

if(length(args) != 6) {stop("Please supply 6 inputs")}

library(data.table)
library(tidyverse)
library(boot)

pheno_id <- args[1]
blocks_file=args[2]
code_file=args[3]
maf_t=as.numeric(args[4])
pval_t=as.numeric(args[5])
output_file=args[6]
print(maf_t)
print(pval_t)
print(pheno_id)

# Download summary stats 
cmd <- paste0("wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv.bgz -O data/raw/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv.bgz")
system(cmd)

# Unzip summary stats 
cmd <- paste0("gunzip -c data/raw/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv.bgz > data/unzipped/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv")
system(cmd)

# Read in summary statistics data
ss_file <- paste0("data/unzipped/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv")
df <- fread(ss_file)

# Check if its a "varorder" file
if (nrow(df) == 0) {
   print("Varorder")
   
   # Download summary stats
    cmd <- paste0("wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/", pheno_id, ".gwas.imputed_v3.both_sexes.varorder.tsv.bgz -O data/raw/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv.bgz")
    system(cmd)

    # Unzip summary stats
    cmd <- paste0("gunzip -c data/raw/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv.bgz > data/unzipped/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv")
    system(cmd)

    # Read in summary statistics data
    ss_file <- paste0("data/unzipped/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv")
    df <- fread(ss_file)
}    
if (nrow(df) == 0) {
    print("Female")

    # Download summary stats
    cmd <- paste0("wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/", pheno_id, ".gwas.imputed_v3.female.tsv.bgz -O data/raw/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv.bgz")
    system(cmd)

    # Unzip summary stats
    cmd <- paste0("gunzip -c data/raw/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv.bgz > data/unzipped/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv")
    system(cmd)

    # Read in summary statistics data
    ss_file <- paste0("data/unzipped/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv")
    df <- fread(ss_file)
}
if (nrow(df) == 0) {
   print("Male")

   # Download summary stats
    cmd <- paste0("wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/", pheno_id, ".gwas.imputed_v3.male.tsv.bgz -O data/raw/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv.bgz")
    system(cmd)

    # Unzip summary stats
    cmd <- paste0("gunzip -c data/raw/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv.bgz > data/unzipped/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv")
    system(cmd)

    # Read in summary statistics data
    ss_file <- paste0("data/unzipped/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv")
    df <- fread(ss_file)

}



# Return the name of the phenotype
file_split <- str_split(ss_file, "\\.")
name_split <- str_split(file_split[[1]][1], "/")
phenotype_code <- name_split[[1]][length(name_split[[1]])]

# Read in the code to phenotype conversion file
all_phenotype_codes <- fread(code_file, header=F)
colnames(all_phenotype_codes) <- c("code", "phenotype")

# Create table of phenotype names that match the code
pheno_name_table <- all_phenotype_codes %>% filter(code == phenotype_code)
phenotype_name <- as.character(pheno_name_table[1,"phenotype"])

# Remove low confidence variants
df_filtered <- df %>% filter(low_confidence_variant == FALSE)

# Remove non-significant variants
df_sig <- df_filtered %>% filter(pval < pval_t)

# Separated column 1 into chromosome, base pair, reference allele, and alt. allele
df_split <- df_sig %>% separate(variant, c("chr", "bp", "ref", "alt"), sep = ":")

# Filtered out low minor allele frequency variants
df_af_filter <- df_split %>% filter(minor_AF > maf_t) %>% filter(chr!= "X")

if(nrow(df_af_filter) == 0) {
  df_out <- as.data.frame(matrix(c(NA, NA, NA, NA, NA), nrow = 1))
  colnames(df_out) <- c("mean_iaf", "lower_ci", "upper_ci", "SNP_number", "phenotype_name")
  print("No significant SNPs")
} else { 


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

  # minSNP <- nrow(df_af_filter) #Delete this for final implementation
  #df_af_filter <- df_af_filter[1:min(minSNP, 100),] #Delete this line for final implementation
  
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
  
  # Record number of SNPs used to calculate average
  nSNPs <- length(df_inc_freq$inc_AF)
  
  # Find average of the increasing allele frequency
  avg_inc_AF <- sum((df_inc_freq$inc_AF * abs(df_inc_freq$beta))) / sum(abs(df_inc_freq$beta))

  if(nSNPs == 1) {
    df_out <- as.data.frame(matrix(c(avg_inc_AF, avg_inc_AF, avg_inc_AF, 1, phenotype_name), nrow = 1))
    colnames(df_out) <- c("mean_iaf", "lower_ci", "upper_ci", "SNP_number", "phenotype_name")
    print("1 significant SNP")
  } else {
  
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
    df_out <- as.data.frame(matrix(c(avg_inc_AF, lb, ub, nSNPs, phenotype_name), nrow = 1))
    colnames(df_out) <- c("mean_iaf", "lower_ci", "upper_ci", "SNP_number", "phenotype_name")
  }
}

# Save output
fwrite(df_out, output_file, row.names = F, col.names = T, sep = "\t")

# Delete raw and zipped files 
cmd <- paste0("rm data/raw/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv.bgz")
system(cmd)

cmd <- paste0("rm data/unzipped/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv")
system(cmd)


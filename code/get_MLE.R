# Feb 2023

# This script will assign each variant to an ld block and filter to only one 
# significant SNP per block.


args=commandArgs(TRUE)

if(length(args) != 7) {stop("Please supply 7 inputs")} #####

library(data.table)
library(tidyverse)
library(boot)

pheno_id <- args[1]
blocks_file=args[2]
code_file=args[3]
maf_t=as.numeric(args[4])
pval_t=as.numeric(args[5])
output_file=args[6]

phenotype_info_file=args[7]############ to access prev and n_cases/controls

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
  df_out <- as.data.frame(matrix(c(NA, NA, NA, 0, phenotype_name, phenotype_code), nrow = 1))
  colnames(df_out) <- c("MLE", "lower_ci", "upper_ci", "SNP_number", "phenotype_name", "phenotype_code") ## <---
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
  
  # Find sign of minor allele
  df_minor_sign <- df_inc_freq %>% mutate(minor_sign = case_when((alt == minor_allele) ~ beta/abs(beta), (alt != minor_allele) ~ -1*beta/abs(beta)))
  df <- df_minor_sign
  
  # Record number of SNPs used to calculate average
  nSNPs <- length(df$inc_AF)
  
  
##------------------------------------------------------------------------------
  
  ## Calculate Non-centrality parameter 
  calc_NCP <- function(p, Nd, Nc, B, prev) {
    
    # Nd = number of cases
    # Nc = number of controls 
    # B = absolute value of liability scale effect size
    # prev = prevalence 
    # p = risk increasing allele freq 
    
    # Compute threshold 
    t <- qnorm(1 - prev)
    
    # Compute mean liability 
    mu_1 <- B * (1 - p)
    mu_2 <- -B * p
    sigma_2 <- 1 - (B^2 * p *(1 - p))
    
    # Conditional Probabilities 
    p_case_a1 <- 1 - pnorm(t, mean = mu_1, sd = sqrt(sigma_2))
    p_case_a2 <- 1 - pnorm(t, mean = mu_2, sd = sqrt(sigma_2))
    p_control_a1 <- pnorm(t, mean = mu_1, sd = sqrt(sigma_2))
    p_control_a2 <- pnorm(t, mean = mu_2, sd = sqrt(sigma_2))
    
    # Compute allele frequency conditional on status
    p_d1 <- (p_case_a1 * p) / (1 - pnorm(t, mean = 0, sd = 1))
    p_c1 <- (p_control_a1 * p) / pnorm(t, mean = 0, sd = 1)
    
    # Observed counts
    q1_o <- 2 * Nd * p_d1
    q2_o <- 2 * Nc * p_c1
    q3_o <- 2 * Nd * (1 - p_d1)
    q4_o <- 2 * Nc * (1 - p_c1)
    
    # Expected counts
    q1_e <- ((2 * (Nd*p_d1 + Nc*p_c1)) * (2*Nd)) / (2* (Nd + Nc))
    q2_e <- ((2 * (Nd*p_d1 + Nc*p_c1)) * (2*Nc)) / (2* (Nd + Nc))
    q3_e <- ((2 * (Nd*(1-p_d1) + Nc*(1- p_c1))) * (2*Nd)) / (2* (Nd + Nc))
    q4_e <- ((2 * (Nd*(1-p_d1) + Nc*(1- p_c1))) * (2*Nc)) / (2* (Nd + Nc))
    
    # Compute NCP 
    NCP <- ((q1_o - q1_e)^2 / q1_e) + ((q2_o - q2_e)^2 / q2_e) + ((q3_o - q3_e)^2 / q3_e) + ((q4_o - q4_e)^2 / q4_e)
    
    # Return NCP
    return(NCP)
  }
  
  # Function to calculate the probability that the minor allele is trait increasing, given its frequency and the selection coefficient
  prob_sign_evo <- function(x, gamma) {
    
    # x = minor allele frequency
    # gamma = population scale selection coefficient
    
    result <- 1 / (1 + exp(gamma * (2*x - 1)))
    return(result)
  }
  
  prob_discovery <- function(x, sign, n_cases, n_controls, B, alpha, prev) {
    
    # x = minor allele frequency
    # sign = sign of liability scale effect size
    # B = absolute value of liability scale effect size
    # alpha =  type 1 error (5e-8 for GWAS)
    # n_cases = number of cases 
    # n_controls = number of controls 
    # prev = prevalence 
    
    # If flip to appropriate increasing allele frequency for given sign 
    if (sign == 1) {
      p <- x
    } else if (sign == -1)  {
      p <- 1 - x 
    } else {
      stop("Please enter a valid sign")
    }
    
    # Calculate ncp 
    ncp <- calc_NCP(p=p, Nd=n_cases, Nc=n_controls, B=B, prev=prev)
    
    # Calculate power 
    power <- 1 - pchisq(qchisq(1- alpha, df =1),df=1, ncp = ncp)
    #print("printing power")
    #print(power)
    
    return(power)
  }
  
  compute_joint_likelihood <- function(gamma, vec_x, obs_sign, B, n_cases, n_controls, alpha, prev) {
    
    # vec_x = vector of minor allele frequencies
    # obs_sign = observed sign
    
    # n_cases = number of cases 
    # n_controls = number of controls
    # B = vector of absolute value of liability scale effect size
    # alpha =  type 1 error (5e-8 for GWAS)
    # prev = prevalence
    # gamma = population scaled selection coefficient
    
    L <- length(vec_x)
    vec_like <- rep(0, L)
    for(i in 1:L) {
      
      P1 <- prob_sign_evo(vec_x[i], gamma) * prob_discovery(x=vec_x[i], sign=1,n_cases=n_cases, n_controls = n_controls, B=B[i], alpha = alpha, prev = prev)
      P2  <- (1 - prob_sign_evo(vec_x[i], gamma)) * prob_discovery(x=vec_x[i], sign=-1,n_cases=n_cases, n_controls = n_controls, B=B[i], alpha = alpha, prev = prev)
      
      if (obs_sign[i] == 1) {
        vec_like[i] <- P1 / (P1 + P2)
      } else if (obs_sign[i] == -1) {
        vec_like[i] <- P2 / (P1 + P2)
      } else {
        stop("Enter valid signs")
      }
    }
    
    final_lik <- sum(log(vec_like))
    return(final_lik)
    
  }
  
##------------------------------------------------------------------------------  
  
  ## Return to the vector of results depending on number of SNPs
  
  phenos_original <- fread(phenotype_info_file)
  phenos <- filter(phenos_original, phenotype == phenotype_code)
  n_cases <- phenos$n_cases
  n_controls <- phenos$n_controls
  prev <- n_cases / (n_cases + n_controls)
  
  if(nSNPs == 1) {
    
    ## Find MLE w/o boot strapping
    
    find_MLE <- function(vec_x, obs_sign, n_cases, n_controls, B, alpha, prev){
      
      MLE <- optimize(f = compute_joint_likelihood, c(0, 10), tol = 0.0001, vec_x = vec_x, obs_sign = obs_sign, n_cases = n_cases, n_controls = n_controls, B = B, alpha = alpha, prev = prev, maximum = TRUE)
      MLE <- MLE$maximum
      
      return(MLE)
    } 
    
    MLE <- find_MLE(vec_x = df$minor_AF, obs_sign = df$minor_sign, B = abs(df$beta), ## Does this need to be absolute
                    n_cases = n_cases, n_controls = n_controls, alpha = pval_t, prev = prev)
    
    df_out <- as.data.frame(matrix(c(MLE, MLE, MLE, 1, phenotype_name, phenotype_code), nrow = 1))
    colnames(df_out) <- c("MLE", "lower_ci", "upper_ci", "SNP_number", "phenotype_name", "phenotype_code")
    print("1 significant SNP")
    
  } else {
    
    ## Boot strapping Method
    MLE_calc <- function(obs_sign, i, vec_x, n_cases, n_controls, prev, B, alpha) {
      
      vec_x <- vec_x[i]
      obs_sign <- obs_sign[i]
      B <- B[i]
      
      MLE <- optimize(f = compute_joint_likelihood, c(0, 10), tol = 0.0001, vec_x = vec_x, obs_sign = obs_sign, n_cases = n_cases, n_controls = n_controls, B = B, alpha = alpha, prev = prev, maximum = TRUE)
      MLE <- MLE$maximum
      return(MLE)
    }
    
    find_MLE_boot <- function(vec_x, obs_sign, n_cases, n_controls, B, alpha, prev){
      
      MLE <- optimize(f = compute_joint_likelihood, c(0, 10), tol = 0.0001, vec_x = vec_x, obs_sign = obs_sign, n_cases = n_cases, n_controls = n_controls, B = B, alpha = alpha, prev = prev, maximum = TRUE)
      MLE <- MLE$maximum
      
      boot_obj <- boot(obs_sign, MLE_calc, R = 1000, n_cases = n_cases, n_controls = n_controls, prev = prev, B = B, alpha = alpha, vec_x = vec_x)
      ci <- boot.ci(boot_obj, type = "norm") ### Figure out what setting to use
      if (is.null(ci)) {
        lc <- NA
        uc <- NA
      } else {
        lc <- ci$normal[2]
        uc <- ci$normal[3]   
      }
      return(c(MLE, lc, uc))
    }
    
    results <- find_MLE_boot(vec_x = df$minor_AF, obs_sign = df$minor_sign, B = abs(df$beta), n_cases = n_cases, n_controls = n_controls, alpha = pval_t, prev = prev)
   
    #### Not getting a confidence interval
    
    
    # Create output data for trait
    MLE <- results[1]
    lb <- results[2]
    ub <- results[3]
    df_out <- as.data.frame(matrix(c(MLE, lb, ub, nSNPs, phenotype_name, phenotype_code), nrow = 1))
    colnames(df_out) <- c("MLE", "lower_ci", "upper_ci", "SNP_number", "phenotype_name", "phenotype_code")
  }
}

##------------------------------------------------------------------------------

# Save output
fwrite(df_out, output_file, row.names = F, col.names = T, sep = "\t")

# Delete raw and zipped files 
cmd <- paste0("rm data/raw/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv.bgz")
system(cmd)

cmd <- paste0("rm data/unzipped/", pheno_id, ".gwas.imputed_v3.both_sexes.tsv")
system(cmd)


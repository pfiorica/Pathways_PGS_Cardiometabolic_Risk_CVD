# Running PRSMix(+)
# Peter Fiorica 
# 15 January 2025

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(PRSmix))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(purrr))

"%&%"=function(a,b) paste(a,b,sep="")

parser <- ArgumentParser(description = "Run combine_PRS with specified arguments")

parser$add_argument("--pheno_file", type = "character", help = "Path to the phenotype file", required = TRUE)
parser$add_argument("--covariate_file", type = "character", help = "Path to the covariate file", required = TRUE)
parser$add_argument("--pgs_file", type = "character", help = "Path to the PGS score file", required = TRUE)
parser$add_argument("--pgs_list", type = "character", help = "Path to the trait-specific PGS list", required = TRUE)
parser$add_argument("--pheno_name", type = "character", help = "Name of the phenotype", required = TRUE)
parser$add_argument("--out", type = "character", help = "Path for output files", required = TRUE)
#parser$add_argument("--covar_list", type = "character", nargs = "+", help = "List of covariates", required = TRUE)
#parser$add_argument("--cat_covar_list", type = "character", nargs = "+", help = "List of categorical covariates", required = TRUE)
parser$add_argument("--isbinary", type = "logical", help = "Whether the phenotype is binary (TRUE/FALSE)", default = TRUE)
#parser$add_argument("--liabilityR2", type = "logical", help = "Compute liability R2 (TRUE/FALSE)", default = TRUE)
parser$add_argument("--IID_pheno", type = "character", help = "Identifier for phenotype column", default = "IID")
#parser$add_argument("--pval_thres_list", type = "numeric", help = "List of p-value thresholds", default = 0.05/435)
#parser$add_argument("--is_extract_adjSNPeff", type = "logical", help = "Extract adjusted SNP effects (TRUE/FALSE)", default = TRUE)
parser$add_argument("--original_beta_files_list", type = "character", help = "The vector contains directories to SNP effect sizes used to compute original PRSs (as weight_file argument from compute PRS above", default = TRUE)
#parser$add_argument("--pop_similarity_file", type = "character", help = "Corresponds to the output from PGSC_CALC ancestry analysis", default = FALSE)
parser$add_argument("--treatment", type = "character", help = "What is the covariate to remove", default = "NA")

args <- parser$parse_args()

treatment <- args$treatment

#################
#COVARIATE
categorical_vars <- c("bl_inc_cat2", "bl_educlvl_5cat", "bl_meno_status", "smoke_status_6m", "prev_anycmd","prev_anycvd", "rad_tx_yn", "anthra_yn", "tras_yn", "horm_yn")

categorical_vars<- setdiff(categorical_vars, treatment)
# Continuous variables (e.g., age, BMI, counts)
continuous_vars <- c( "dxage", "bmi_new",  "pc_util_count")

outcomes <- c("inc_anycmd", "inc_hyperten", "inc_diabetes", "inc_dyslipid", "inc_anycvd", "inc_seriouscvd", "inc_arrhythmia", "inc_hf", "inc_ischemic", "inc_stroke", "inc_vte", "inc_cvddeath")
cat_def <- c("M","M", "M","M", "V","V","V","V", "V", "V","V", "V")
dict <-setNames(cat_def, outcomes)

if(dict[[args$pheno_name]]=="M"){
  categorical_vars <- categorical_vars[categorical_vars !="prev_anycmd"]
}else{
  categorical_vars <-  categorical_vars[categorical_vars !="prev_anycvd"]
}

covariates <- c(categorical_vars, continuous_vars)

scores_file_vector <-fread(args$pgs_file, header =F)$V1 

pgs_list <- fread(args$pgs_list, header = F)
num_pgs <- nrow(pgs_list)

print("Processing the phenotype: " %&% args$pheno_name %&%" ; Number of PGS Corresponding to it: " %&% num_pgs)

if(num_pgs > 1) {
  print("The number of PGS for this phenotype is greater than one, so we are proceeding with `combined_PRS()`.")
  combine_PRS(
    pheno_file = args$pheno_file ,
    covariate_file = args$covariate_file ,
    score_files_list = scores_file_vector ,
    trait_specific_score_file = args$pgs_list ,
    pheno_name = args$pheno_name ,
    isbinary = TRUE , 
    out = args$out, 
    liabilityR2 = TRUE ,
    IID_pheno =  args$IID_pheno ,
    covar_list = covariates,
    cat_covar_list = categorical_vars,
    is_extract_adjSNPeff = FALSE,
    train_size_list = NULL,
    #power_thres_list = 0.95,
    pval_thres_list = (0.05/7),
    read_pred_training = FALSE,
    read_pred_testing = FALSE,
    original_beta_files_list = args$original_beta_files_list
  )
}else{
  print("The number of PGS for this phenotype is only one, so we are proceeding with `eval_single_PRS()`.")
  phenotype <- fread(args$pheno_file, header = T)
  covariate_file<-fread(args$covariate_file, header = T)
  pgs_name <- paste0(pgs_list$V1[1],"_SUM", sep="")
  PRS <- fread(scores_file_vector[1], header = T) %>% select(FID, IID, pgs_name)
  df_to_join <- list(phenotype, covariate_file, PRS)
  
  combined_file <- reduce(df_to_join,left_join, by = c("FID", "IID"))
  combined_file <- combined_file %>% 
    select(args$pheno_name, pgs_name, all_of(covariates))
  pgs_data <- c()
  single_pgs <- eval_single_PRS(
    data_df= combined_file,
    pheno = args$pheno_name ,
    prs_name = pgs_name,
    covar_list = covariates,
    isbinary = TRUE,
    liabilityR2 = TRUE,
    alpha = 0.05
  )
  
  fwrite(single_pgs, args$out %&% "_power.0.95_pthres.0.00714285714285714_test_summary_trait_eval_single.txt" )
}


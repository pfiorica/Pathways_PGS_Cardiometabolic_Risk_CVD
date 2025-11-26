# PXSTools Script
# 24 January 2025
# Peter Fiorica

# Libraries
library(tidyverse)
library(data.table)
library(ggsci)
library(PXStools)
library(rlog)
library(pROC)
library(boot)

#### FUNCTIONS
# Concatenate Function
"%&%"=function(a,b) paste(a,b,sep="")

# Bootstrapping function
boot_auc <- function(data, indices) {
  seed = 1
  d <- data[indices, ]
  roc(d$phenotype, d$pred)$auc
}

# Define input files

covariate_file <- fread("/bd-fs-mnt/TenantShare/exthome/home/l578311/pgs_sandbox/phenotypes/ALL_covariates.txt", header = T)
pheno_list <- fread("/bd-fs-mnt/TenantShare/exthome/home/l578311/sig_phenotypes.txt", header = F)$V1
#pheno_list<- c("inc_anycmd", "inc_hyperten", "inc_diabetes", "inc_dyslipid", "inc_anycvd", "inc_seriouscvd", "inc_arrhythmia", "inc_hf", "inc_ischemic", "inc_stroke", "inc_vte", "inc_cvddeath")
phenotypes <- fread("/bd-fs-mnt/TenantShare/exthome/home/l578311/pgs_sandbox/phenotypes/large_phenotype_dummy.txt", header = T)
pgs_list <- fread("/bd-fs-mnt/TenantShare/exthome/home/l578311/janise_share/pgs_list.txt", header = F)$V1
pgs_file <- fread("/bd-fs-mnt/TenantShare/exthome/home/l578311/janise_share/pw_32_pgs_Zscore.txt", header = T)
out_path <- "/bd-fs-mnt/TenantShare/exthome/home/l578311/pgs_sandbox/"
pgs_list_dir <- "/bd-fs-mnt/TenantShare/exthome/home/l578311/prsmix_lists/"
fdr_threshold <- 0.05

# Make output directory if it does not exist
directory_path <- out_path

# Check if the directory exists, and create it if it doesn't
if (!dir.exists(directory_path %&% "/pxs_output")) {
  dir.create(directory_path %&% "/pxs_output", recursive = TRUE)
  message("Directory created: ", directory_path %&% "/pxs_output")
} else {
  message("Directory already exists: ", directory_path %&% "/pxs_output")
}
covariate_file <- covariate_file %>% dplyr::select(-ancestry)

# Initial Manipulation
covariates <- colnames(covariate_file)[3:length(colnames(covariate_file))]
combined_file <- left_join(phenotypes, pgs_file, by = c("ID"="IID"))
pgs_cols <- paste0(pgs_list, "_SUM")

pgs_cols <- gsub( "_hmPOS_GRCh37", "", pgs_cols)

characterization <- c("CMD", "CMD", "CMD" , "CMD", "CVD", "CVD", "CVD", "CVD")
pheno_list <- data.table(pheno_list, type=characterization)

# Split Ratios
group_a_rat <- 0.6 # First ratio for XWAS
group_b_rat <- 0.2
group_c_rat <- 0.2

## Calculate sizes
n <- length(covariate_file$IID)
train_size <- floor(group_a_rat * n)
test_size <- floor(group_b_rat * n)
validation_size <- n - train_size - test_size
set.seed(123)

#Functional Aspect
print("Processing the Exposure-wide Association Study")
summary_res <- data.table()
summary_res_pgs <- data.table()

for(pheno in pheno_list$pheno_list){
  print("Beginning to process " %&% pheno %&% " for XWAS...")
  characterization <- pheno_list %>%
    dplyr::filter(pheno_list == pheno) %>%
    dplyr::pull(type)
  
  ######## SPLITTING INDIVIDUALS TO THREE GROUPS
  # Shuffle IDs
  #shuffled_IDs <- sample(covariate_file$IID)
  # Assign splits
  #group_a_IDs <- shuffled_IDs[1:train_size]
  #group_b_IDs <- shuffled_IDs[(train_size + 1):(train_size + test_size)]
  #group_c_IDs <- shuffled_IDs[(train_size + test_size + 1):n]
  
  XWAS_res <- data.table()
  # SELECT APPROPRIATE VARIABLES
  df <- combined_file %>% select(ID,
                                 PHENO=!!sym(pheno),
                                 all_of(covariates)) %>%
    mutate(PHENO=if_else(is.na(PHENO), 0, PHENO)) %>%
    filter(PHENO %in% c("0","1"))
  
  #Assure that unknowns are represented
  unknown_rows <- which(df$bl_educlvl_5cat=="Unknown" | df$smoke_status_6m =="Unknown")
  unknown_df <- df[unknown_rows,]
  rest_df<-df[-unknown_rows,]
  
  ## Calculate sizes
  n <- nrow(rest_df)
  unique_IDs<- unique(rest_df$ID)
  
  # Define Training Data
  train_ids <- sample(unique_IDs, size = group_a_rat* length(unique_IDs))
  remain_ids <- setdiff(unique_IDs, train_ids)
  
  # Define validation data
  val_ids <- sample(remain_ids, size = 0.5*length(remain_ids))
  test_ids <- setdiff(remain_ids, val_ids)
  ######## SPLITTING INDIVIDUALS TO THREE GROUPS
  # Assign the unknowns

  unknown_ids <- unique(unknown_df$ID)
  shuffled_ids <- sample(unknown_ids)
  n <- length(shuffled_ids)
  n_train <- floor(n / 3)
  n_val   <- floor((n - n_train) / 2)
  n_test  <- n - n_train - n_val
  train_unknown_ids <- shuffled_ids[1:n_train]
  val_unknown_ids   <- shuffled_ids[(n_train + 1):(n_train + n_val)]

  test_unknown_ids  <- shuffled_ids[(n_train + n_val + 1):n]
  
  # Assign splits
  group_a_IDs <- c(train_ids, train_unknown_ids)
  group_b_IDs <- c(test_ids, test_unknown_ids)
  group_c_IDs <- c(val_ids, val_unknown_ids)
  
  if(characterization=="CMD"){
    COV <- c("dxage", "prev_anycvd")
    XVAR <- setdiff(covariates, COV)
    REM="L"
    print("Performing XWAS for phenotpe: " %&% pheno)
  }else{
    COV <- c("dxage", "prev_anycmd")
    XVAR <- setdiff(covariates, COV)
    REM="L"
    print("Performing XWAS for phenotpe: " %&% pheno)
  }
  ######## PERFORM XWAS
  XWAS_res <- xwas(df=df,
                   X= XVAR,
                   cov=COV,
                   mod = "logistic",
                   IDA = group_a_IDs,
                   removes = REM)
  fwrite(XWAS_res, out_path %&% "/pxs_output/" %&% pheno %&% "_pxstools_" %&% length(XVAR)  %&%"_variables_out.txt",  col.names = T, row.names = T, sep = "\t", quote = F)
  
  #PGS XWAS 
  df_pgs <- combined_file %>% select(ID,
                                     PHENO=!!sym(pheno),
                                     all_of(covariates),
                                     pgs_cols) %>%
    mutate(PHENO=if_else(is.na(PHENO), 0, PHENO)) %>%
    filter(PHENO %in% c("0","1"))
  select_pgs<-fread(pgs_list_dir %&% pheno %&% "_PGS_list.txt", header = F)$V1
  select_pgs <- paste0(select_pgs, "_SUM")
  if(characterization=="CMD"){
    COV <- c("dxage", "prev_anycvd")
    XVAR_pgs <- setdiff(covariates, COV)
    XVAR_pgs <- c(XVAR_pgs, select_pgs)
    REM="L"
    print("Performing XWAS for phenotpe: " %&% pheno)
  }else{
    COV <- c("dxage", "prev_anycmd")
    XVAR_pgs <- setdiff(covariates, COV)
    XVAR_pgs <- c(XVAR_pgs, select_pgs)
    REM="L"
    print("Performing XWAS for phenotpe: " %&% pheno)
  }
  
  ######## PERFORM XWAS
  XWAS_pgs <- xwas(df=df_pgs,
                   X= XVAR_pgs,
                   cov=COV,
                   mod = "logistic",
                   IDA = group_a_IDs,
                   removes = REM)
  fwrite(XWAS_res, out_path %&% "/pxs_output/" %&% pheno %&% "_pxstools_" %&% length(XVAR)  %&%"_variables_out_pgs.txt",  col.names = T, row.names = T, sep = "\t", quote = F)
  
  ### 
  #POLY-EXPOSURE SCORE
  print("Processing the Poly-exposure_score_analysis...")
  exposures_to_keep <- XWAS_res %>% filter(fdr < fdr_threshold)
  exposures_rownames <- rownames(exposures_to_keep)
  XVAR_sig <- XVAR[sapply(XVAR, function(x) any(grepl(x, exposures_rownames)))]
  
  # Run poly-exposure Score
  print("Calculating Poly-exposure score... NOT PGS")
  PXSS = PXSgl(df=df,
               X= XVAR_sig,
               cov=COV,
               mod = "logistic",
               IDA = group_a_IDs,
               IDB = group_b_IDs,
               IDC=group_c_IDs,
               removes = REM,
               seed = 1)
  
  
  ### PERFORMING FOR pgs
  exposures_to_keep_pgs <- XWAS_pgs %>% filter(fdr < fdr_threshold)
  exposures_rownames_pgs <- rownames(exposures_to_keep_pgs)
  XVAR_sig_pgs <- XVAR_pgs[sapply(XVAR_pgs, function(x) any(grepl(x, exposures_rownames_pgs)))]
  
  # Run poly-exposure Score
  print("Calculating Poly-exposure score... FOR PGS")
  PXS_pgs = PXSgl(df=df_pgs,
                  X= XVAR_sig_pgs,
                  cov=COV,
                  mod = "logistic",
                  IDA = group_a_IDs,
                  IDB = group_b_IDs,
                  IDC=group_c_IDs,
                  removes = REM,
                  seed = 1)
  
  
  # Bootstrap CI for AUC  
  print("Calculating AUCs for PXS ...")
  
  roc_obj <- roc(PXSS$PHENO, PXSS$pred)
  auc_value <- auc(roc_obj)
  print("PXS AUC of ROC " %&% pheno %&% ": " %&% auc_value)
  
  data <- data.frame(pred = PXSS$pred, phenotype = PXSS$PHENO)
  results <- boot(data, boot_auc, R = 1000)
  ci <- boot.ci(results, type = "bca")
  lower_ci <- ci$bca[4]
  upper_ci <- ci$bca[5]
  
  
  print("Calculating AUCs for PXS ...for PGS")
  
  roc_obj_pgs <- roc(PXS_pgs$PHENO, PXS_pgs$pred)
  auc_value_pgs <- auc(roc_obj_pgs)
  print("PXS AUC of ROC " %&% pheno %&% ": " %&% auc_value_pgs)
  
  data_pgs <- data.frame(pred = PXS_pgs$pred, phenotype = PXS_pgs$PHENO)
  results_pgs <- boot(data_pgs, boot_auc, R = 1000)
  ci_pgs <- boot.ci(results_pgs, type = "bca")
  lower_ci_pgs <- ci_pgs$bca[4]
  upper_ci_pgs <- ci_pgs$bca[5]
  
  # Save plot
  png(out_path %&% "/pxs_output/" %&% pheno %&% "_roc_curve.png")
  plot(
    roc_obj,
    col = "blue",
    main = "ROC Curve Comparison for " %&% pheno,
    xlim = c(1, 0),  # Set x-axis range
    ylim = c(0, 1),  # Set y-axis range
    xaxs = "i",      # Remove padding from x-axis
    yaxs = "i",      # Remove padding from y-axis
    lwd = 2          # Set line width
  )
  # Plot PGS Data
  plot(
    roc_obj_pgs,
    col = "red",
    lty = 2,         # Dashed line for PGS
    lwd = 2,
    add = TRUE       # Add to the existing plot
  )
  
  # Add hard axes at 0-1
  abline(h = 0, col = "black")  # Horizontal line at y=0
  abline(h = 1, col = "black")  # Horizontal line at y=1
  abline(v = 0, col = "black")  # Vertical line at x=0
  abline(v = 1, col = "black")  # Vertical line at x=1
  
  # Annotate the AUC on the plot
  # Annotate the AUC values
  auc_label <- paste("Non-PGS AUC =", round(auc_value, 3))  # Format AUC value for non-PGS
  auc_label_pgs <- paste("PGS AUC =", round(auc_value_pgs, 3))  # Format AUC value for PGS
  text(0.5, 0.15, auc_label, col = "blue", cex = 1.2)  # Position and style for non-PGS AUC
  text(0.5, 0.1, auc_label_pgs, col = "red", cex = 1.2) 
  
  legend(
    "bottomright",
    legend = c("Non-PGS", "PGS"),
    col = c("blue", "red"),
    lty = c(1, 2),   # Solid line for non-PGS, dashed for PGS
    lwd = 2,
    cex = 1.2        # Increase text size for readability
  )
  dev.off()
  
  # Store results in summary data.table
  summary_res <- rbind(summary_res, data.table(
    phenotype = pheno,
    auc = auc_value,
    ci_lower = lower_ci,
    ci_upper = upper_ci,
    XVAR_sig = paste(XVAR_sig, collapse = ", ")
  ))
  summary_res_pgs <- rbind(summary_res_pgs, data.table(
    phenotype = pheno,
    auc = auc_value_pgs,
    ci_lower = lower_ci_pgs,
    ci_upper = upper_ci_pgs,
    XVAR_sig = paste(XVAR_sig_pgs, collapse = ", ")
  ))
}



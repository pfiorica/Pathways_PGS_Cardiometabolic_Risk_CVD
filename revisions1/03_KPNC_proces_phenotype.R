# Process Phenotype_files

library(dplyr)
library(tidyr)
library(data.table)
library(readxl)

"%&%"=function(a,b) paste(a,b,sep="")

data.dir <-"/bd-fs-mnt/TenantShare/exthome/home/l578311/cvd/pgs_revised/data/final_cohorts/"

data_files <-list.files(path=data.dir, full.names = T)
geno_map <- read_excel("/bd-fs-mnt/TenantShare/exthome/home/l578311/cvd/cvd/PW_genotype_ID_092823.xlsx")
geno_map <-geno_map %>%
  mutate(STUDYID= gsub(" 3268", "03268", STUDYID))%>%
  mutate(STUDYID= sprintf("%05s", STUDYID))
rds_list <- lapply(data_files, readRDS)

rds_list <- lapply(rds_list, as.data.table)

names(rds_list) <- tools::file_path_sans_ext(basename(data_files))


overall <- rds_list[["pgs_zscores_overall_peter"]]

pheno_start<-overall


pheno_no_pgs <- pheno_start %>% select(-starts_with("z_pgs_"))

pheno <- pheno_no_pgs  %>%
  rename(STUDYID = studyid)
  
pheno <- pheno %>%
    left_join(geno_map, by = c("STUDYID")) %>% select(-scanID)


# Edit BMI



#pheno <- update_BMI(pheno)

# Edit for INCOME:
update_column_with_dictionary <- function(data, column_name, dictionary) {
  # Ensure the column exists
  if (!column_name %in% colnames(data)) {
    stop("Column not found in the dataset: ", column_name)
  }
  # Convert the column to character (if not already)
  data[[column_name]] <- as.character(data[[column_name]])
  # Map values using the dictionary
  data[[column_name]] <- sapply(data[[column_name]], function(x) {
    if (x %in% names(dictionary)) {
      dictionary[[x]]
    } else {
      x # Keep original value if no match in the dictionary
    }
  })
  return(data)
}

income_mapping <- list(
  "25K-49K" = "II",
  "50K-89K" = "III",
  "<25K" = "I",
  ">=90K" = "IV",
  "Unknown" = "X"
)

pheno<- update_column_with_dictionary(pheno, "bl_inc_cat2", income_mapping)
pheno$bl_educlvl_5cat <- gsub(" ", "_",pheno$bl_educlvl_5cat)
pheno$bl_educlvl_5cat <- gsub("t-g", "tG",pheno$bl_educlvl_5cat)


outcomes <- c("inc_anycmd", "inc_hyperten", "inc_diabetes", "inc_dyslipid", "inc_anycvd", "inc_seriouscvd", "inc_arrhythmia", "inc_hf", "inc_ischemic", "inc_stroke", "inc_vte", "inc_cvddeath")
cat_def <- c("M","M", "M","M", "V","V","V","V", "V", "V","V", "V")
covariates<- c("dxage", "bl_inc_cat2", "bl_educlvl_5cat", "bmi_new", "bl_meno_status", "smoke_status_6m", "pc_util_count", "prev_anycmd", "prev_anycvd", "rad_tx_yn", "anthra_yn", "tras_yn", "horm_yn")

dict <-setNames( cat_def, outcomes)


pheno <- pheno 
dir<- "/bd-fs-mnt/TenantShare/exthome/home/l578311/pgs_sandbox/phenotypes/"
#Write one big phenotype file
fwrite(pheno, dir %&% "large_phenotype_dummy.txt", col.names = T, row.names = F, sep = "\t", quote = F)

# Write individual phenotype files
for (i in outcomes){
  a <- pheno %>%
    mutate(IID=ID, FID=ID) %>%
    select(IID, FID, !!sym(i)) %>%
    mutate(!!sym(i) := if (all(is.na(!!sym(i)) | !!sym(i) == 1)) {
      replace_na(!!sym(i), 0)
    } else {
      !!sym(i)
    })
  a <- a %>% filter(!!sym(i) %in% c("0","1"))
  b <- pheno %>% 
    mutate(IID=ID, FID=ID) %>%
    select(IID, FID, all_of(covariates), ancestry)  %>%
    mutate(bmi_new=ifelse(IID=="02994","28.52", bmi_new)) %>%
    filter(IID %in% a$IID) %>%
    filter(!is.na(bmi_new))
  ref<-dict[[i]]
  if(ref=="M"){
    b<-b %>% select(-prev_anycmd)
  }else{
    if (ref=="V") {
      b<-b %>% select(-prev_anycvd)
    }
  }
  fwrite(b, dir %&% "ALL_" %&% i %&% "_covariates.txt", col.names = T, row.names = F, sep = "\t", quote = F )
  fwrite(a, dir %&% "ALL_" %&% i %&% "_phenotype.txt",  col.names = T, row.names = F, sep = "\t", quote = F)
}


# Write one big covariate file

big_covar <- pheno %>% 
  mutate(IID=ID, FID=ID) %>%
  select(IID, FID, all_of(covariates), ancestry)  %>%
  mutate(bmi_new=ifelse(IID=="02994","28.52", bmi_new)) %>%
  filter(!is.na(bmi_new))

fwrite(big_covar, dir %&% "ALL_covariates.txt", col.names = T, row.names = F, sep = "\t", quote = F)


### Ancestry Specific Files


for (i in outcomes){
  for (anc in unique(pheno$ancestry)){
    a <- pheno %>%
      filter(ancestry == anc) %>%
      mutate(IID=ID, FID=ID) %>%
      select(IID , FID , !!sym(i)) %>%
      mutate(!!sym(i) := if (all(is.na(!!sym(i)) | !!sym(i) == 1)) {
        replace_na(!!sym(i), 0)
      } else {
        !!sym(i)
      })
    a <- a %>% filter(!!sym(i) %in% c("0","1"))
    b <- pheno %>%
      mutate(IID=ID, FID=ID) %>%
      select(IID, FID, all_of(covariates)) %>%
      filter(IID %in% a$IID) %>%
      filter(!is.na(bmi_new))
    ref<-dict[[i]]
    if(ref=="M"){
      b<-b %>% select(-prev_anycmd)
    }else{
      b<-b %>% select(-prev_anycvd)
    }
    fwrite(b, dir %&% anc %&% "_" %&% i %&% "_covariates.txt", col.names = T, row.names = F, sep = "\t", quote = F )
    fwrite(a, dir %&% anc %&% "_" %&% i %&%"_phenotype.txt",  col.names = T, row.names = F, sep = "\t", quote = F)
  }
}

#for (anc in unique(pheno$ancestry)){
#  b <- pheno %>%
#    filter(ancestry == anc) %>%
#    mutate(IID=ID, FID=ID) %>%
#    select(IID, FID, all_of(covariates))
#  fwrite(b, dir %&% anc %&% "_" %&% "_covariates.txt", col.names = T, row.names = F, sep = "\t", quote = F )
#}

# Lastly, we want to write the PGS lists...

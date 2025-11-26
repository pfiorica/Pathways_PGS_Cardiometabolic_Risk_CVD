# Pull List for PRSMix+
# Peter Fiorica
#17 January 2025

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(grid))

"%&%"=function(a,b) paste(a,b,sep="")

format_string <- function(x) {
  x <- tolower(x)        # Convert to lowercase
  x <- gsub(" ", "_", x) # Replace spaces with underscores
  return(x)
}

setwd("/bd-fs-mnt/TenantShare/exthome/home/l578311/")

sheet_name <- "/bd-fs-mnt/TenantShare/exthome/home/l578311/cvd/pgs_revised/data/results_data/all_results_norounding.rds"

total <- readRDS(sheet_name)


PGS_list<-fread("janise_share/Selected_PGS_Labeled.csv", header = T)
PGS_list$Trait<-format_string(PGS_list$Trait)
PGS_list$Trait<-paste0("z_pgs_",PGS_list$Trait)
PGS_labeled<-PGS_list %>% select(Trait, Category, `PGS #`)

total<- total %>%
  rename(Trait=term)#, p_white = p_eur, p_black=p_afr, p_hispanic=p_amr, p_asian=p_eas)

working_names <- c("inc_anycmd", "inc_hyperten", "inc_diabetes", "inc_dyslipid", "inc_anycvd", "inc_seriouscvd", "inc_arrhythmia", "inc_hf", "inc_ischemic", "inc_stroke", "inc_vte", "inc_cvddeath")
formal_names<- c("Any CMD risk factor", "Hypertension","Diabetes","Dyslipidemia", "Any CVD",  "Serious CVD", "Arrhythmia","Heart failure or cardiomyopathy","Ischemic heart disease","Stroke","VTE","CVD-related death")


map_key<-setNames(unlist(working_names), unlist(formal_names))
#outcome_rename<-fread("reference_files/survival_outcome_rename.csv", header = T)



total <- total %>%
  #left_join(outcome_rename, by = c("Survival_Outcome" = "Old_Names")) %>%
  #mutate(Survival_Outcome = coalesce(Working_Names, Survival_Outcome)) %>%
  #select(-Working_Names) %>% 
  filter(!Trait %in% c("z_pgs_cad2", "z_pgs_cad3", "z_pgs_cardiovascular disease", "z_pgs_systolic_blood_pressure", "z_pgs_type_2_diabetes")) %>%
  mutate(outcome=map_key[as.character(outcome)])

PGS_total_labeled <- left_join(total, PGS_labeled, by = "Trait") %>%
  mutate(`PGS #`=if_else(Trait == "z_pgs_cad3725", "PGS003725",
                         if_else(Trait =="z_pgs_hb_a1c","PGS000685",
                                 if_else(Trait =="z_pgs_diabetes", "PGS002308",
                                         if_else(Trait == "z_pgs_sbp", "PGS002807", `PGS #`)))))

#PGS_total_labeled <- total %>% mutate(`PGS #`=Trait)



output_dir <- "prsmix_lists"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Calculate the Bonferroni correction threshold
bonferroni_threshold <- 0.05 / 324 #nrow(PGS_total_labeled)

# Filter associations with P.value_Overall < Bonferroni threshold
significant_associations <- PGS_total_labeled %>%
  filter(p_overall < bonferroni_threshold) %>%
  mutate(Survival_Outcome = gsub(" ", "_", outcome))

# Create an empty list to store file paths
file_paths_list <- list()

# For each survival outcome, get the list of significant `PGS #` and generate file paths
survival_outcomes <- unique(significant_associations$Survival_Outcome)

for (outcome1 in survival_outcomes) {
  # Filter the significant associations for the specific survival outcome
  outcome_data <- significant_associations %>%
    filter(Survival_Outcome == outcome1)
  
  # Extract the unique `PGS #` for this outcome
  pgs_numbers <- unique(outcome_data$`PGS #`)
  
  # Generate file paths for each `PGS #`
  file_paths <- pgs_numbers %>%
    paste0("/bd-fs-mnt/TenantShare/exthome/home/l578311/janise_share/pgs_results/", 
           ., "_zscore_score_pgsc_calc.txt")
  
  # Store the file paths in the list
  file_paths_list[[outcome1]] <- file_paths
  
  # Write the `PGS #` list to a text file
  writeLines(pgs_numbers, file.path(output_dir, paste0(outcome1, "_PGS_list.txt")))
  
  # Write the file paths to a text file
  writeLines(file_paths, file.path(output_dir, paste0(outcome1, "_file_paths.txt")))
}

# Save the file paths list as an RDS object (optional for later use)
saveRDS(file_paths_list, file.path(output_dir, "file_paths_list.rds"))

# Output to indicate script completion
message("Script completed. Significant associations and file paths written to '", output_dir, "'.")


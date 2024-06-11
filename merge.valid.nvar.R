#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# idx <- as.numeric(args[1])
# paste0("idx=",idx)
# paste0("class of idx is ",class(idx))

rm(list = ls())
dir_PGS <- "/projects/rpci/songyao/Haiyang/PGS/SBP.Diabetes.2PGS_08122023/weights"
# dir_PGS <- "C:/Users/Haiyang/Box/RA/RoswellPark/2022/1010 CMD and CVD analysis in PWS/PGS/weights"
PGS.files = list.files(dir_PGS, full.names = T)
# PGS.scoreID <- substr(PGS.files, start = 59, stop = 67)
PGS.scoreID <- unlist(regmatches(PGS.files, gregexpr("(PGS00.*?)(?=\\.txt)", PGS.files, perl = TRUE)))
# PGS.scoreNM <- gsub(" ", ".", substr(PGS.files, start = 69, stop = nchar(PGS.files)-4))
PGS.scoreNM <- c("Diabetes", "SBP")

df_PGS.sm_Pathway <- as.data.frame(matrix(NA, ncol = 8, nrow = length(PGS.scoreNM)))
names(df_PGS.sm_Pathway) <- c("Trait", "N.of.Variants", "N.of.Ambg", "% Ambg", "Total.Matched", 
                            "% Total.Matched", "Flipped.Mathced", "% Flipped.Matched")

# df.sm_ls<-list()
for (i in 1:length(PGS.scoreID)) {
  for (pop in c("AA")) {
  # for (pop in c("AA", "EA4", "EU", "HIS")) {
    mywd = paste("/projects/rpci/songyao/Haiyang/PGS/SBP.Diabetes.2PGS_08122023/", PGS.scoreID[i], "/",  pop, "/", sep = '')
    # mywd = paste("C:/Users/Haiyang/Box/RA/RoswellPark/2022/1010 CMD and CVD analysis in PWS/PGS/", PGS.scoreID[i], "/",  pop, "/", sep = '')
    setwd(mywd)
    cat(getwd(), sep = '\n')
    
    # trait
    df_PGS.sm_Pathway$Trait[i] <- PGS.scoreNM[i]
    
    # n of variants
    weight <- read.table(PGS.files[i], header = T, sep = "\t", as.is = T, 
                         na.strings = c(" ", "", "NA", NA), fill = T)
    df_PGS.sm_Pathway$N.of.Variants[i] <- nrow(weight)
    
    # n of ambg and %
    if ("other_allele" %in% names(weight)) {names(weight)[names(weight) == "other_allele"] <- "reference_allele"}
    df_PGS.sm_Pathway$N.of.Ambg[i] <- sum(paste0(weight$effect_allele, ":", weight$reference_allele) %in% c("A:T", "T:A", "C:G", "G:C"))
    df_PGS.sm_Pathway$`% Ambg`[i] <- round(100 * df_PGS.sm_Pathway$N.of.Ambg[i] / nrow(weight), 2)
    
    # total matched and %
    nvar.files <- list.files()
    nvar.files <- nvar.files[grepl("nvar", nvar.files)]
    nvar.files <- nvar.files[grepl("valid", nvar.files)]
    nvar.list <- lapply(nvar.files, function(x) {read.table(x, header = F, sep = "\t")})
    # for (j in 1:length(score.list)) {
    #   names(score.list[[j]]) <- c("samples", substr(score.files[j], 1, nchar(score.files[j])-10))
    # }
    nvar <- sum(do.call("rbind", nvar.list))
    df_PGS.sm_Pathway$Total.Matched[i] <- nvar
    df_PGS.sm_Pathway$`% Total.Matched`[i] <- round(100*nvar/nrow(weight), 2)
    # flipped matched and %
    flip.files <- list.files()
    flip.files <- flip.files[grepl("unmatched", flip.files)]
    flip.list <- lapply(flip.files, function(x) {read.table(x, header = F, sep = "\t")})
    nflip <- sum(do.call("rbind", flip.list))
    df_PGS.sm_Pathway$Flipped.Mathced[i] <- nflip
    df_PGS.sm_Pathway$`% Flipped.Matched`[i] <- round(100*nflip/nrow(weight), 2)
    
    # all.score$total<-rowSums(all.score[,-which(names(all.score)=="samples")])
    # assign(paste0("all.score.", pop), all.score)
  }
  # df_ls[[i]]<-rbind(all.score.AA[,c('samples', "total")], all.score.EA4[,c('samples', "total")], 
  #                   all.score.EU[,c('samples', "total")], all.score.HIS[,c('samples', "total")])
  # names(df_ls[[i]])[2]<-paste0("PGS_",PGS.scoreNM[i])
  cat(paste0("i = ", i, "; ", PGS.scoreNM[i], " is done."), sep = '\n')
}

# df_PGS_Pathway<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "samples", all.x = TRUE),
#                        df_ls)

library(openxlsx)
write.xlsx(df_PGS.sm_Pathway, "/projects/rpci/songyao/Haiyang/PGS/SBP.Diabetes.2PGS_08122023/PGS.matching.summary_08142023.xlsx", colWidths = "auto")
# write.csv(df_PGS.sm_Pathway, "/projects/rpci/songyao/Haiyang/PGS/60PGS/PGS.matching.summary_01132022.csv", row.names = F)


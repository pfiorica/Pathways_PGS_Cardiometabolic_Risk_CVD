#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# idx <- as.numeric(args[1])
# paste0("idx=",idx)
# paste0("class of idx is ",class(idx))

dir_PGS <- "/projects/rpci/songyao/Haiyang/PGS/SBP.Diabetes.2PGS_08122023/weights"
PGS.files = list.files(dir_PGS, full.names = T)
# PGS.scoreID <- substr(PGS.files, start = 59, stop = 67)
PGS.scoreID <- unlist(regmatches(PGS.files, regexpr("PGS00.*(?=\\.txt)", PGS.files, perl = TRUE))) 
# PGS.scoreNM <- gsub(" ", ".", substr(PGS.files, start = 69, stop = nchar(PGS.files) - 4))
PGS.scoreNM <- c("Diabetes", "SBP")
  
# df_PGS_Pathway<-as.data.frame(matrix(NA, ncol = length(PGS.scoreNM)+1, nrow = 1))
# names(df_PGS_Pathway)<-c("samples", paste0("PGS_",PGS.scoreNM))

df_ls <- list()
for (i in 1:length(PGS.scoreID)) {
  for (pop in c("AA", "EA4", "EU", "HIS")) {
    mywd = paste("/projects/rpci/songyao/Haiyang/PGS/SBP.Diabetes.2PGS_08122023/", PGS.scoreID[i], "/",  pop, "/", sep = '')
    setwd(mywd)
    cat(getwd(), sep = '\n')
  
    score.files <- list.files()
    score.files <- score.files[!grepl("n", score.files)]
    score.list <- lapply(score.files, function(x) {read.table(x, header = F, sep = "\t")})
    for (j in 1:length(score.list)) {
      names(score.list[[j]]) <- c("samples", substr(score.files[j], 1, nchar(score.files[j]) - 10))
    }
    all.score <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "samples", all.x = TRUE),
                      score.list)
    all.score$total <- rowSums(all.score[, -which(names(all.score) == "samples")])
    assign(paste0("all.score.", pop), all.score)
  }
  df_ls[[i]] <- rbind(all.score.AA[, c('samples', "total")], all.score.EA4[, c('samples', "total")], 
                    all.score.EU[, c('samples', "total")], all.score.HIS[, c('samples', "total")])
  names(df_ls[[i]])[2] <- paste0("PGS_", PGS.scoreNM[i])
  cat(paste0("i=", i, "; ", PGS.scoreNM[i], " is done."), sep = '\n')
}

df_PGS_Pathway <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "samples", all.x = TRUE),
                         df_ls)
library(openxlsx)
write.xlsx(df_PGS_Pathway, "/projects/rpci/songyao/Haiyang/PGS/SBP.Diabetes.2PGS_08122023/PGS_08142023.xlsx", colWidths = "auto")
write.csv(df_PGS_Pathway, "/projects/rpci/songyao/Haiyang/PGS/SBP.Diabetes.2PGS_08122023/PGS_08142023.csv", row.names = F)


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
p <- as.numeric(args[1])
paste0("p = ", p)
paste0("class of p is ",class(p))

# dir_PGS <- "C:/Users/Haiyang/Box/RA/RoswellPark/2022/1010 CMD and CVD analysis in PWS/PGS/weights"
dir_PGS <- "/projects/rpci/songyao/Haiyang/PGS/SBP.Diabetes.2PGS_08122023/weights"
PGS.files = list.files(dir_PGS, full.names = T)
# PGS.scoreID <- substr(PGS.files, start = 59, stop = 67)
PGS.scoreID <- unlist(regmatches(PGS.files, regexpr("PGS00.*(?=\\.txt)", PGS.files, perl = TRUE)))
PGS.nms <- c("Diabetes", "SBP")

# Print the result

dir.create(paste("./", PGS.scoreID[p], sep = ''))
dir.create(paste("./", PGS.scoreID[p], "/TotalScores", sep = ''))

# the target GWAS data with weights
# Polygenic score information where the first rows are:
#rsID, Chr, position, effect allele, other allele, effect weight
weightFile = PGS.files[p]
weight = read.table(weightFile, header = T, sep = "", as.is = T, na.strings = c(" ", "", "NA", NA), fill = T)
  if ("other_allele" %in% names(weight)) {names(weight)[names(weight) == "other_allele"] <- "reference_allele"}
  head(weight)
  # in this weight file the effect allele and reference allele are the same with what we have
  # weight$effect_allele_rv<-weight$reference_allele
  # weight$reference_allele_rv<-weight$effect_allele
  weight$name = paste(weight$chr_name, weight$chr_position, weight$effect_allele, weight$reference_allele, sep=":")
  # The last row is NA
  # weight=weight[1:143,]

library(GWASTools)
# library(survival)
# library(cmprsk)

# # chromosome number
# chr = rep(1:22,4)[idx]
# 
# # Race population
# pop = rep(c("AA", "EA4", "EU", "HIS"), each=22)[idx]

for (pop in c("AA", "EA4", "EU", "HIS")) {
  dir.create(paste("./", PGS.scoreID[p], "/", pop, sep = ''))
  cat(paste0("pop ", pop), sep = "\n")
  for (chr in (1:22)[1:22 %in% unique(weight$chr_name)]) {
    cat(paste0("chr", chr), sep = "\n")
    
    # output file for current chromosome and population
    outfile = paste("./", PGS.scoreID[p], "/", pop, "/chr.", chr, "scores.txt", sep = '')
    
    # directory of the imputation data
    gdsDir = paste("/projects/rpci/qzhu/eschille/BreastGWAS/No_Y_Chrom/maf1.noHM/Input/Imputation/", pop,
                   "cohort/Sample.Filtered/1kg/R2filter/GDS", sep = '')
    # gdsDir = "C:/Users/Haiyang/Box/RA/RoswellPark/2021/7.28 PGS Calculation"
    # scan all the file names with "chr "current chrnumber"-.dgs" files
    files = list.files(gdsDir, pattern = paste("chr", chr, "-.*.gds$", sep = ""), full.names = T)
    # output files for number of variants and valid variants
    outfile1 = paste(outfile, "nvar", sep = ".")
    outfile2 = paste(outfile, "valid.nvar", sep = ".")
    outfile3 = paste(outfile, "n.unmatched.alleles", sep = ".")
    
    scores = c()
    nvar = 0
    valid.nvar = 0
    umchd.Alle = 0
    # for each file of current population and chromosome
    for (i in 1:length(files)) {
      gdsfile = files[i]
      snpAnnotFile = sub("gds", "snp.RData", gdsfile) 
      
      gds <- GdsGenotypeReader(gdsfile)
      # add phenotype information
      samples = getVariable(gds, "sample.id") # or samples = getScanID(gds) 
      load(snpAnnotFile)
      # Get the snps and alleles to match with the weights data
      anno = pData(snpAnnot)
      
      genoData <-  GenotypeData(gds)
      # the imputation data, the large matrix to multiply with the weights
      geno <- getGenotype(genoData, use.names = T)
      close(gds)
      
      #variants = apply(anno[,c("chromosome", "position", "alleleA", "alleleB")], 1, function(x){paste(c(as.numeric(x[1:2]), sort(x[3:4])), collapse=":")})
      variants1 = paste(anno$chromosome, anno$position, anno$alleleB, anno$alleleA, sep=":")
      variants2 = paste(anno$chromosome, anno$position, anno$alleleA, anno$alleleB, sep=":")
      ind1 = match(variants1, weight$name)
      ind2= match(variants2, weight$name)
      ind1[is.na(ind1)]<-ind2[is.na(ind1)]
      ind <- ind1
      nv = sum(!is.na(ind))
      if (nv==0) {
        cat(paste("No variant match in", gdsfile, sep=" ") , sep="\n")
        next
      }
      if(nv!=0){
        nvar = nvar + nv 
        t = which(!is.na(ind))
        anno = anno[t,]
        geno = geno[t,]
        weight1 = weight[ind[t],]
        
        valid = nrow(weight1[!is.na(weight1$effect_weight) & weight1$effect_weight != 0,])
        valid.nvar = valid.nvar + valid 
        
        ind = which(anno$alleleB != weight1$effect_allele)
        umchd.Alle = umchd.Alle + length(ind)
        if(length(ind) == 0){
          if(ncol(data.frame(geno)) > 1){
            # scores is N of cases by matched snps
            scores = cbind(scores, t(geno) %*% weight1$effect_weight)}
          if(ncol(data.frame(geno)) == 1){
            scores = cbind(scores, geno * weight1$effect_weight)}
        }
        
        if(length(ind) > 1){
          # imputation 0~2, refallele not matched, use 2-the imputated value to 
          # change the direction
          geno[ind,] = 2 - geno[ind,]
          scores = cbind(scores, t(geno) %*% weight1$effect_weight)
        }
        
        if(length(ind) == 1){
          if(ncol(data.frame(geno)) == 1) {
            geno = 2-geno
            scores = cbind(scores, geno * weight1$effect_weight)
          } else {
            geno[ind,] = 2 - geno[ind,]
            scores = cbind(scores, t(geno) %*% weight1$effect_weight)
          }
        }
        
      }
    }
    
    scores = rowSums(scores, na.rm = T)
    write.table(scores, file = outfile, sep = "\t", col.names = F, row.names = T, quote = F, append = F)
    write.table(nvar, file = outfile1, sep = "\t", col.names = F, row.names = F, quote = F, append = F)
    write.table(valid.nvar, file=outfile2, sep="\t", col.names = F, row.names = F, quote = F, append = F)
    write.table(umchd.Alle, file=outfile3, sep="\t", col.names = F, row.names = F, quote = F, append = F)
    
  }
}




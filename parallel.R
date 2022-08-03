if(FALSE){ #Unused code
args <- commandArgs(trailingOnly = TRUE)

#For developmental accessibilty data
gwas_sum <- readRDS('/external/rprshnas01/kcni/jxia/chromatin-accessibility/GWAS_sumstats/SZ_gwas_final_list.rds')
#To check if code works, create a fake entry
gwas_sum[nrow(gwas_sum)+1,] <- list("fake_developmental", 1, 10090, "B", "C", 1)

Developmental_DAR <- readRDS('/external/rprshnas01/kcni/jxia/chromatin-accessibility/Developmental/Developmental_DAR.rds')

row <- gwas_sum[args,]
chrom <- row[["CHR"]] #Stores the location of each SNP 
loc <- row[["BP"]]
for (i in 1:length(Developmental_DAR)){
  if (Developmental_DAR[[i,"chr"]]==chrom & Developmental_DAR[[i,"start"]]<=loc & Developmental_DAR[[i,"end"]]>=loc){
    row <- cbind(row, Developmental_DAR[i,])
    saveRDS(as.list(row), file = paste("parallel_results/SZ_intersected_",as.character(args),".rds", sep=""))
  }
}
}








args <- commandArgs(trailingOnly = TRUE)

DAR_grouped <- readRDS("~/chromatin-accessibility/Adult/DAR_grouped.rds")
gwas_sum <- readRDS('/external/rprshnas01/kcni/jxia/chromatin-accessibility/GWAS_sumstats/SZ_gwas_final_list.rds')

gwas_sum[nrow(gwas_sum)+1,] <- list("fake_adult", 13, 20231045, "B", "C", 1)
#For adult chromatin accessibility data
row <- as.list(gwas_sum[args,])

chrom <- row[["CHR"]] #Stores the location of each SNP 
loc <- row[["BP"]]
for (j in 1:length(DAR_grouped)){ #Loops through each cell subtype
  subtype <- names(DAR_grouped[j])
  row[subtype] <- 0.00
  for (k in 1:nrow(DAR_grouped[[j]])){ #Loops through each accessible region of the subtype
    if(DAR_grouped[[j]][[k,"chr"]]==chrom & DAR_grouped[[j]][[k,"start"]]<=loc & DAR_grouped[[j]][[k,"end"]]>=loc){
      row[subtype] <- DAR_grouped[[j]][k,"log.fold.change."] #If the SNP is within the accessible region, store the log FC in 
      next
    }
  }
}

saveRDS(row, file = paste("parallel_results/SZ_intersected_",as.character(args),".rds", sep=""))




#qbatch -w 00:10:00 --ppj 1 xx.txt
#ppj is the number of cores
#showq
#cancel jobnumber
args <- commandArgs(trailingOnly = TRUE)

gwas_sum <- readRDS('/external/rprshnas01/kcni/jxia/chromatin-accessibility/BD_gwas_final_list.rds')
DAR_grouped <- readRDS('/external/rprshnas01/kcni/jxia/chromatin-accessibility/DAR_grouped.rds')


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

saveRDS(row, file = paste("parallel_results/BD_intersected_",as.character(args),".rds", sep=""))

#qbatch -w 00:10:00 --ppj 1 xx.txt
#ppj is the number of cores
#showq
#cancel jobnumber
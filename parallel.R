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







#Argument passed in from shell script
args <- as.numeric(commandArgs(trailingOnly = TRUE))
ending <- args+4
#Import chromatin accessibility & gwas data
DAR_grouped <- readRDS("~/chromatin-accessibility/Adult/DAR_grouped.rds")
gwas_sum <- readRDS('/external/rprshnas01/kcni/jxia/chromatin-accessibility/GWAS_sumstats/SZ_gwas_finemapped_list.rds')
#Add fake entry to gwas list to verify that code is working
gwas_sum[nrow(gwas_sum)+1,] <- list("fake_adult", 13, 20231045, "B", "C", 1)
gwas_sum[args:ending,"accessible_region"] <- "unavailable"
df <- data.frame()

for (i in args:ending){
  #Stop looping if beyond the number of gwass
  if (i > nrow(gwas_sum)){
    break
  }
  if (i==args){
    #For adult chromatin accessibility data
    row <- as.list(gwas_sum[i,])
    chrom <- row[["CHR"]] #Stores the location of each SNP 
    loc <- row[["BP"]]
    for (j in 1:length(DAR_grouped)){ #Loops through each cell subtype
      subtype <- names(DAR_grouped[j])
      row[subtype] <- 0.00
      #row["accessible_region"] <- "unavailable"
      for (k in 1:nrow(DAR_grouped[[j]])){ #Loops through each accessible region of the subtype
        if(DAR_grouped[[j]][[k,"chr"]]==chrom & DAR_grouped[[j]][[k,"start"]]<=loc & DAR_grouped[[j]][[k,"end"]]>=loc){
          row[subtype] <- DAR_grouped[[j]][k,"log.fold.change."] #If the SNP is within the accessible region, store the log FC in 
          row["accessible_region"] <- DAR_grouped[[j]][k,"Location"]
          next
        }
      }
    }
    df <- rbind(df,row)
    #df[nrow(df)+1,] <- row
  } else{
    #For adult chromatin accessibility data
    row <- as.list(gwas_sum[i,])
    chrom <- row[["CHR"]] #Stores the location of each SNP 
    loc <- row[["BP"]]
    for (j in 1:length(DAR_grouped)){ #Loops through each cell subtype
      subtype <- names(DAR_grouped[j])
      row[subtype] <- 0.00
      #row["accessible_region"] <- "unavailable"
      for (k in 1:nrow(DAR_grouped[[j]])){ #Loops through each accessible region of the subtype
        if(DAR_grouped[[j]][[k,"chr"]]==chrom & DAR_grouped[[j]][[k,"start"]]<=loc & DAR_grouped[[j]][[k,"end"]]>=loc){
          row[subtype] <- DAR_grouped[[j]][k,"log.fold.change."] #If the SNP is within the accessible region, store the log FC in 
          row["accessible_region"] <- DAR_grouped[[j]][k,"Location"]
          next
        }
      }
    }
    #df <- rbind(df,row)
    df[nrow(df)+1,] <- row
    
  }
  
}
saveRDS(df, file = paste("parallel_results/SZ_intersected_",as.character(args),".rds", sep=""))





#qbatch -w 00:10:00 --ppj 1 xx.txt
#ppj is the number of cores
#showq
#cancel jobnumber
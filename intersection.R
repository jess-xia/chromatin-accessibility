library(dplyr)
library(janitor)
library(readxl)
library(stringr)
library(data.table)

#To check the datatype of each column
#str()
#To restart R session
#.rs.restartR()

#Read in supplementary Table 14c. SNARE-Seq2 Differentially Accessible Regions (DARs, q value < 0.001 and log(fold change) > 1) identified by subclass for human M1
#Specify n_max= to read in only specified number of rows
DAR <- read_xlsx("~/DAR.xlsx") %>%
  rename_with(make.names) %>%
  arrange(., desc(log.fold.change.)) %>%
  #Extract chromosome & start/end location and insert as a separate column
  mutate(chr=str_match(Location, "chr([A-Z0-9]{1,2}):")[,2], 
         start=as.numeric(str_match(Location, ":([0-9]+)-([0-9]+)")[,2]),
         end=as.numeric(str_match(Location, ":([0-9]+)-([0-9]+)")[,3]),
         .after=Location)
DAR_grouped <- split(DAR, f=DAR$Subclass)


#Loads scz gwas sumstats. SNP=RSid, A1=effect allele, A2=other allele, beta=effect estimate, P=p-value
SZ <- readRDS("/external/rprshnas01/kcni/jxia/chromatin-accessibiliy/SZ_gwas_final_list.rds") 
BD <- readRDS("/external/rprshnas01/kcni/jxia/chromatin-accessibiliy/BD_gwas_final_list.rds")

#For each GWAS SNP, loop through the differential accessible regions by cell subtype. If there is an overlap, the log fold change is store as extra columns appended to the original gwas sumstat file
for (i in 17:nrow(BD)){ #Loops through each SNP
  chrom <- BD[[i, "CHR"]] #Stores the location of each SNP 
  loc <- BD[[i, "BP"]]
  for (j in 1:length(DAR_grouped)){ #Loops through each cell subtype
    subtype <- names(DAR_grouped[j])
    BD[i,subtype] <- 0.00
    for (k in 1:nrow(DAR_grouped[[j]])){ #Loops through each accessible region of the subtype
      if(DAR_grouped[[j]][[k,"chr"]]==chrom & DAR_grouped[[j]][[k,"start"]]<=loc & DAR_grouped[[j]][[k,"end"]]>=loc){
        BD[i,subtype] <- DAR_grouped[[j]][k,"log.fold.change."] #If the SNP is within the accessible region, store the log FC in 
        next
      }
    }
  }
}

#BD <- subset(BD, select = c("SNP","CHR", "BP", "A1","A2","P"))


'''
Function to perform intersection. Problem is nothing is saved when stopping prematurely. Better to not use the function?
intersection <- function(gwas_sum){
  for (i in 1:nrow(gwas_sum)){ #Loops through each SNP
    chrom <- gwas_sum[[i, "CHR"]] #Stores the location of each SNP 
    loc <- gwas_sum[[i, "BP"]]
    for (j in 1:length(DAR_grouped)){ #Loops through each cell subtype
      subtype <- names(DAR_grouped[j])
      gwas_sum[i,subtype] <- 0.00
      for (k in 1:nrow(DAR_grouped[[j]])){ #Loops through each accessible region of the subtype
        if(DAR_grouped[[j]][[k,"chr"]]==chrom & DAR_grouped[[j]][[k,"start"]]<=loc & DAR_grouped[[j]][[k,"end"]]>=loc){
          gwas_sum[i,subtype] <- DAR_grouped[[j]][k,"log.fold.change."] #If the SNP is within the accessible region, store the log FC in 
          next
        }
      }
    }
    return(gwas_sum)
  }
}
intersected_BD <- intersection(BD)
'''
  
  
#Adding a fake DAR
#scz_withloc[nrow(scz_withloc)+1,] <- list("fakesnp", NA, NA, NA, NA, NA, 13, 20231046, NA, NA)
#Removing the fake DAR. The negative value in head() specifies how many rows to remove from the bottom of the dataframe
#scz_withloc <- head(scz_withloc, -1)
#Delete the SNP column in the DAR data table
#DAR <- subset(DAR, select=-SNP)
  
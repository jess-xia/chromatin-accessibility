

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
  mutate(chr=as.numeric(str_match(Location, "chr([0-9]{1,2}):")[,2]), 
         start=as.numeric(str_match(Location, ":([0-9]+)-([0-9]+)")[,2]),
         end=as.numeric(str_match(Location, ":([0-9]+)-([0-9]+)")[,3]),
         .after=Location)
DAR_grouped <- split(DAR, f=DAR$Subclass)


#Reference file that matches RSids to chromosomal locations for GWAS summary statistics
rsid_ref <- fread("/external/rprshnas01/netdata_kcni/stlab/Alex_PRS/annotation/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz",
                  select = c("#CHROM", "POS", "ID", "REF", "ALT"))
rsid_ref <- rename(rsid_ref, "SNP"="ID") #Rename rsid column to be the same as in the summary statistics. Get error related to not enough memory when combined with above


#Load BD gwas sumstats if chromosomal location not matched to rsid yet
#bd <- fread("~/BD.sumstats") %>%
#  left_join(., rsid_ref, by="SNP")
#Load BD gwas sumstats if chromosomal location already matched to rsid 
bd <- fread("~/BD.sumstats")
#fwrite(bd, file = "~/BD_intersected.sumstats",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


#For each GWAS SNP, loop through the differential accessible regions by cell subtype. If there is an overlap, the log fold change is store as extra columns appended to the original gwas sumstat file
for (i in 1:nrow(bd)){ #Loops through each SNP
  chrom <- bd[i,`#CHROM`] #Stores the location of each SNP 
  loc <- bd[i, POS]
  for (j in 1:length(DAR_grouped)){ #Loops through each cell subtype
    subtype <- names(DAR_grouped[j])
    bd[i,subtype] <- 0.00
    for (k in 1:nrow(DAR_grouped[[j]])){ #Loops through each accessible region of the subtype
      if(DAR_grouped[[j]][[k,"chr"]]==chrom & DAR_grouped[[j]][[k,"start"]]<=loc & DAR_grouped[[j]][[k,"end"]]>=loc){
        bd[i,subtype] <- DAR_grouped[[j]][k,"log.fold.change."] #If the SNP is within the accessible region, store the log FC in 
        next
      }
    }
  }
}

#Loads scz gwas sumstats. SNP=RSid, A1=effect allele, A2=other allele, beta=effect estimate, P=p-value
scz <- read.table("C:/Tripathy Lab/Intersection/SZ.sumstats") %>%
  row_to_names(row_number = 1) %>% #Reassigns first row as column names using janitor package
  transform(., P = as.numeric(P)) %>% #Converts p values from chr to numeric type
  mutate(., log = -log10(P)) %>% #Add a column with -log10(P)
  filter(., log >= 8) %>% #Keep only rows with -log10(P)
  arrange(., desc(log)) %>% #Arrange SNPs by decreasing -log10(P) value
  left_join(., rsid_ref, by="SNP")

#Intersects gwas summary with rsid file, assigns chromosomal location for each rsid entry in gwas summary statistics. fread is for regular delimited files; i.e., where every row has the same number of columns. delimited flat file contains one or more records set off from each other by a specified delimiter, or separator.
scz_withloc <- fread("/external/rprshnas01/netdata_kcni/stlab/Alex_PRS/annotation/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz",
                     select = c("#CHROM", "POS", "ID", "REF", "ALT")) %>%
  rename("SNP"="ID") %>% #Rename rsid column to be the same as in scz
  
  
  #Adding a fake DAR
  #scz_withloc[nrow(scz_withloc)+1,] <- list("fakesnp", NA, NA, NA, NA, NA, 13, 20231046, NA, NA)
  #Removing the fake DAR. The negative value in head() specifies how many rows to remove from the bottom of the dataframe
  #scz_withloc <- head(scz_withloc, -1)
  #Delete the SNP column in the DAR data table
  #DAR <- subset(DAR, select=-SNP)
  
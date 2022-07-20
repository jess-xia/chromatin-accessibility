

library(data.table)
library(dplyr)
library(h2o)

#Sumstat preparation
#NOTE: Different GWAS papers provides the raw sumstat file in different formats
setwd('/external/rprshnas01/netdata_kcni/stlab/Alex_PRS/sumstats/raw_sumstats') #Original files are stored here
standard_colnames=c("SNP","CHR","A1","A2","P") #Standardized column names to be used for each GWAS


#Schizophrenia (original paper:'https://doi.org/10.1038/s41586-022-04434-5')
SZ <- fread("PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz") %>%
  mutate(., freq=(FCAS*NCAS+FCON*NCON)/(NCAS+NCON)) %>% #Calculate average frequency of the risk allele across cases and controls
  mutate(., MAF=ifelse(freq > .5, 1-freq, freq)) %>% #If freq of allele of interest indicates its the major allele, subtract from 1 to calculate minor allele frequency (MAF), otherwise leave as is. 
  filter(., MAF>=0.01 & PVAL<5*10^-8 & IMPINFO>=0.8) %>% #Filter out SNPs with MAF<1% (standard)
  select(., c('ID','CHROM','A1','A2','PVAL')) #Columns kept in final dataframe

colnames(SZ)=standard_colnames #To standardize the names of columns
fwrite(SZ, file = "YOURFOLDER/SZ.sumstats",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


#Bipolar disorder (original paper:'https://doi.org/10.1038/s41588-021-00857-4')
BD <- fread("pgc.bip.full.2012-04.txt.gz") #Error when loading, addressed in next line of code
colnames(BD) <- c("SNP","CHR","BP","A1","A2","OR","SE","P","INFO","NGT","CEUaf") #Fixes the error

ref <- fread("reference.1000G.maf.0.005.txt.gz",
             select = c("SNP", "MAF")) #CEUaf=frequency of A1 but is missing many values. Use this reference instead

BD <- BD %>% left_join(.,ref,by="SNP") %>% #Combine ref to BD sumstats raw file
  filter(., MAF>=0.01 & P<5*10^-8 & INFO>=0.8) %>%
  select(.,c("SNP","CHR","A1","A2","P"))

fwrite(BD, file = "~/BD.sumstats",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)



#Depression (original paper:'https://doi.org/10.7488/ds/2458')
MDD<-fread("PGC_UKB_depression_genome-wide.txt") %>%
  mutate(., MAF=ifelse(Freq > .5, 1-Freq, Freq)) %>% #Convert Freq to MAF, same as before
  filter(., MAF>=0.01 & P<5*10^-8) %>% #There's no INFO column this time, only filtering by MAF & P value
  select(., c("MarkerName","A1","A2","P"))

MDD$A1 <- toupper(MDD$A1) #Convert A1 and A2 to uppercase like the other dataframes
MDD$A2 <- toupper(MDD$A2)

colnames(MDD)=standard_colnames
fwrite(MDD, file = "YOURFOLDER/MDD.sumstats",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)



######with these you should be able to prepare IQ, EA, CP and OCD
IQ<-fread("SavageJansen_2018_intelligence_metaanalysis.txt")
#doi='https://doi.org/10.1038/s41588-018-0152-6'

EA<-fread("EA4_additive_excl_23andMe.txt.gz")
#doi='https://doi.org/10.1038/s41588-018-0147-3'

CP<-fread("GWAS_CP_all.txt")
#doi='https://doi.org/10.1038/s41588-018-0147-3'

OCD=fread("/external/rprshnas01/netdata_kcni/stlab/Alex_PRS/sumstats/raw_sumstats/ocd_aug2017.gz")
#doi='https://doi.org/10.1038/mp.2017.154'
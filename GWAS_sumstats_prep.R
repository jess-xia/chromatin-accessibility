library(data.table)
library(dplyr)
library(h2o)
library(readxl)
library(stringr)
library(janitor)
library(rtracklayer)

#GWAS summary statistics preparation, standardizes column names for each GWAS, ensures RSID is matched to chromosomal location
#NOTE: Different GWAS papers provides the raw sumstat file in different formats
#setwd('/external/rprshnas01/kcni/jxia/chromatin-accessibiliy') #Original files are stored here
standard_colnames=c("SNP","CHR", "BP", "A1","A2","P") #Standardized column names to be used for each GWAS


#Schizophrenia (original paper:'https://doi.org/10.1038/s41586-022-04434-5')
SZ <- read_xlsx("/external/rprshnas01/kcni/jxia/chromatin-accessibility/GWAS_sumstats/SZ_gwas_sum.xlsx", sheet=2) %>%
  #as.data.frame() %>%
  filter(., `P-comb`<5*10^-8) %>% #P-comb = p-value of combined discover-replication meta analysis. Filter out SNPs with P-value less than 5*10^-8, yields 342 SNPs
  mutate(A1=str_match(A1A2, "([A-Z])/([A-Z])")[,2],
         A2=str_match(A1A2, "([A-Z])/([A-Z])")[,3], 
         liftover=paste("chr", as.character(CHR), ":", as.character(BP), "-", as.character(BP), sep="")) #Formated coordinates for liftover

#Export list of coordinates in the correct format for http://genome.ucsc.edu/cgi-bin/hgLiftOver
fwrite(as.list(SZ$liftover), file="GWAS_sumstats/SZ_gwas_hg19_coordinates.txt", sep="\n")
#New corresponding coordinates
new_coor <- fread("~/chromatin-accessibility/GWAS_sumstats/hglft_genome_1df00_932fa0.bed", header=FALSE)
#Remove SNPs where conversion failed, code varies based on failed conversions
hg38_SZ <- filter(SZ, CHR!=23) %>%
  cbind(., new_coor) %>%
  mutate(CHR=str_match(V1, "chr([A-Z0-9]{1,2}):")[,2], 
         BP=as.numeric(str_match(V1, ":([0-9]+)-([0-9]+)")[,2])) %>% #Extract chromosomal location from new coordinates
  subset(select=c('SNP','CHR', 'BP', 'A1', 'A2','P-comb')) #Columns kept in final dataframe

colnames(hg38_SZ) <- standard_colnames #To standardize the names of columns
saveRDS(hg38_SZ, file = "GWAS_sumstats/SZ_gwas_final_list.rds") #Save as a single R object


#Bipolar disorder (original paper:'https://doi.org/10.1038/s41588-021-00857-4')
BD <- read_xlsx("/external/rprshnas01/kcni/jxia/chromatin-accessibility/BD_gwas_sum.xlsx", sheet=3) %>% #p-values are already filtered for P<5×10^−8
  row_to_names(row_number = 1) %>%
  mutate(A1=str_match(`A1/A2`, "([A-Z])/([A-Z])")[,2],
         A2=str_match(`A1/A2`, "([A-Z])/([A-Z])")[,3]) %>%
  select(., c('SNP','CHR', 'BP', 'A1', 'A2','P')) %>%
  head(.,-5) # Last 5 rows are descriptions

BD[ , c('BP', 'P')] <- apply(BD[ , c('BP', 'P')],2, #2 specifies apply is applied by column
                                    function(x){as.numeric(as.character(x))}) #Convert specified columns to numeric type
colnames(BD) <- standard_colnames
saveRDS(BD, file = "BD_gwas_final_list.rds")


#Depression (original paper:'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6522363/')
MDD <- read_xlsx("/external/rprshnas01/kcni/jxia/chromatin-accessibility/MDD_gwas_sum.xlsx") %>%
  subset(select=c(2,1,3,4,5,31)) %>%
  row_to_names(row_number = 2) %>%
  head(-1) #Last row contains description

colnames(MDD) <- standard_colnames
MDD$A1 <- toupper(MDD$A1) #Convert A1 and A2 to uppercase like the other dataframes
MDD$A2 <- toupper(MDD$A2)
MDD$BP <- as.numeric(as.character(MDD$BP))
MDD$P <- as.numeric(as.character(MDD$P))

saveRDS(MDD, file = "MDD_gwas_final_list.rds")


'''
#Intersects gwas summary with rsid file, assigns chromosomal location for each rsid entry in gwas summary statistics. fread is for regular delimited files; i.e., where every row has the same number of columns. delimited flat file contains one or more records set off from each other by a specified delimiter, or separator.
scz_withloc <- fread("/external/rprshnas01/netdata_kcni/stlab/Alex_PRS/annotation/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz",
                     select = c("#CHROM", "POS", "ID", "REF", "ALT")) %>%
  rename("SNP"="ID") %>% #Rename rsid column to be the same as in scz
#Reference file that matches RSids to chromosomal locations for GWAS summary statistics
rsid_ref <- fread("/external/rprshnas01/netdata_kcni/stlab/Alex_PRS/annotation/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz",
                  select = c("#CHROM", "POS", "ID", "REF", "ALT"))
rsid_ref <- rename(rsid_ref, "SNP"="ID") #Rename rsid column to be the same as in the summary statistics. Get error related to not enough memory when combined with above
'''

######with these you should be able to prepare IQ, EA, CP and OCD
IQ<-fread("SavageJansen_2018_intelligence_metaanalysis.txt")
#doi='https://doi.org/10.1038/s41588-018-0152-6'

EA<-fread("EA4_additive_excl_23andMe.txt.gz")
#doi='https://doi.org/10.1038/s41588-018-0147-3'

CP<-fread("GWAS_CP_all.txt")
#doi='https://doi.org/10.1038/s41588-018-0147-3'

OCD=fread("/external/rprshnas01/netdata_kcni/stlab/Alex_PRS/sumstats/raw_sumstats/ocd_aug2017.gz")
#doi='https://doi.org/10.1038/mp.2017.154'
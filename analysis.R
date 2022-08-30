library(biomaRt)
library(readxl)

SZ_intersected_hg38 <- readRDS("~/chromatin-accessibility/Adult/SZ_intersected_hg38_with_.rds") %>%
  filter(SNP!="fake_adult")
rsids <- unique(SZ_intersected_hg38$SNP)

#Compare gwas hits before and after genome reference liftover
SZ_intersected_mismatched_genome_reference <- read_csv("~/chromatin-accessibility/Adult/SZ_intersected.csv")
rsids_original <- unique(SZ_intersected_mismatched_genome_reference$SNP)
intersect(rsids, rsids_original)

#Import finemapped results of SNP hits
finemapped_SZ <- read_excel("~/chromatin-accessibility/GWAS_sumstats/SZ_finemapped.xlsx", sheet = 2) %>%
  filter(rsid %in% rsids)

#Match associated gene and finemapped results to each rsid
joined <- full_join(SZ_intersected_hg38, finemapped_SZ, by = c("SNP" = "rsid")) %>%
  subset(., select=c("SNP","CHR", "BP", "A1","A2","P", "cell_type", "accessibility", "index_snp", "finemap_posterior_probability", "gene_symbol", "gene_ensembl", "ensembl_gene_classification", "impact"))


#Check if SNP hits are actually in accessible regions
snp_chr = "11"
snp_pos = 134426490
DAR <- read_excel("~/chromatin-accessibility/DAR.xlsx") %>%
  rename_with(make.names) %>%
  arrange(., desc(log.fold.change.)) %>%
  #Extract chromosome & start/end location and insert as a separate column
  mutate(chr=str_match(Location, "chr([A-Z0-9]{1,2}):")[,2], 
         start=as.numeric(str_match(Location, ":([0-9]+)-([0-9]+)")[,2]),
         end=as.numeric(str_match(Location, ":([0-9]+)-([0-9]+)")[,3]),
         .after=Location) %>%
  filter(chr==snp_chr,  start < snp_pos, end > snp_pos)


#Find genes not retrieved from finemapped data
#listMarts()
ensembl_snp <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
#listDatasets(ensembl107)
#View(listAttributes(ensembl_snp))
matched <- getBM(attributes=c("associated_gene", "ensembl_gene_name", 'refsnp_id'), 
      filters = 'snp_filter', 
      values = rsids,
      mart=ensembl_snp) %>%
  subset(select = c("ensembl_gene_name",  "refsnp_id"))

joined <- left_join(joined, matched, by=c("SNP"="refsnp_id"))


#Paper identified 342 linkage-disequilibrium-independent significant SNPs located in 287 loci. Match SNP to loci based on provided table
SZ_342SNP_287loci <- read_excel("~/chromatin-accessibility/GWAS_sumstats/SZ_342SNP_287loci.xls")

for (i in 1:nrow(SZ_intersected_hg38)){
  for (j in 1:nrow(SZ_342SNP_287loci)){
    if(grepl(SZ_intersected_hg38[[i, "SNP"]], SZ_342SNP_287loci[[j, "indices"]], fixed=TRUE)){
      SZ_intersected_hg38[[i, "genes_all"]] <- SZ_342SNP_287loci[[j, "genes_all"]]
      next
    }
  }
}
genes_all <- SZ_intersected_hg38$genes_all

genes_all <- unique(unlist(strsplit(genes_all,",", fixed=TRUE)))



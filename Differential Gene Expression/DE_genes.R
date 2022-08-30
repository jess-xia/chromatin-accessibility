DE_genes_merged_tables_final <- readRDS("~/chromatin-accessibility/Differential Gene Expression/DE_genes_merged_tables_final.rds")
SZ_intersected_hg38 <- read_csv("~/chromatin-accessibility/Adult/SZ_intersected_hg38_with_nearest_gene.csv") %>%
  filter(SNP!="fake_adult")
SZ_intersected_hg38 <- subset(SZ_intersected_hg38, select=-X1)
  

#Convert nested list to dataframe
gene_list <- SZ_intersected_hg38$nearest_gene %>%
  append(c("Nr1H2", "HINFP", "SMAD2", "HOXB3", "NR2C1", "POU5F1", "NR2F2", "Nr1H4", "NKX2-4"))

#gene_list <- genes_all

DE_df <- bind_rows(DE_genes_merged_tables_final, .id = "cell_type") %>%
  filter(gene %in% gene_list) %>%
  filter(Meta_adj.P.Val < 0.01) %>%
  subset(., select=colnames(.)[1:7])

library(data.table)
library(dplyr)
library(readr)
library(tidyverse)
library(janitor)

#Read in files
output <- read.table(file = 'PsyOPS_output.tsv', sep = '\t', header = TRUE)
SZ <- read.csv(file="SZ_intersection_summary.csv")
gene_expression <- read_csv("new_CgG_avg_count.csv")

#Joins Psyops results with previously intersected SNPs
psyops_combined <- left_join(SZ, output, by=c("SNP"="lead_variant"))  %>%
  mutate(mapped_cell_types = case_when(!differentially_accessible_cell_type %in% c("L2.3.IT", "L5.IT", "L6.IT") ~ differentially_accessible_cell_type, 
                                       differentially_accessible_cell_type=="L2.3.IT" ~ "IT", 
                                       differentially_accessible_cell_type=="L5.IT" ~ "IT",
                                       differentially_accessible_cell_type=="L6.IT" ~ "IT"))%>%
  filter(P<5*10^-7 & finemap_posterior_probability>0.1)

#Add additional columns to gene expression data
gene_expression$cell_type_with_highest_expression <- colnames(gene_expression)[max.col(subset(gene_expression, select = -X1), ties.method = "first")+1]
gene_expression <- gene_expression %>%
  clean_names() %>%
  mutate(highest_expression = do.call(pmax, subset(., select = oligodendrocyte:l4_it)))

#IT in gene expression sheet counts for 3 cell types in SNARE-seq dataset: L2.3.IT, L5.IT, L6.IT
#Joins gene expression with previously intersected SNPs
combined <- left_join(psyops_combined, gene_expression, by=c("gene"="x1"))%>%
  rename(Micro.PVM=microglia,
         Oligo=oligodendrocyte, 
         VIP=vip,
         OPC=opc, 
         Astro=astrocyte, 
         PVALB=pvalb, 
         LAMP5=lamp5,
         SST=sst, 
         IT=it) %>%
  pivot_longer(Oligo:l4_it,names_to = "expression_cell_type", values_to = "expression") %>%
  filter(mapped_cell_types==expression_cell_type & expression>1) %>%
  arrange(desc(finemap_posterior_probability), desc(PsyOPS_score), desc(expression)) 


write_csv(subset(combined, select=-X), "combined_notfilteredbypsyops.csv")

#Saves the list of genes in a line-separated text file
cat(unique(combined$gene), sep="\n", file="./gene_list.txt")
cat(colnames(combined), sep="\n")



#Shreejoy's code
snpcells_df = read_csv('combined_notfilteredbypsyops.csv')

marker_df = read_csv('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/All_hodge_regions/Ranked_markers_ALLReg_ITexpand_WL35IT_lfct11_minpct15_dup.csv')

# below, i'm adding a new cell column called new cell name that is just the first 5 characters, this helps me match cell type names
# between jessica's file and sonny's marker data frame
df = left_join(snpcells_df %>% mutate(new_cell_name = str_sub(differentially_accessible_cell_type, start = 1, end = 5), 
                                      new_cell_name = case_when(differentially_accessible_cell_type == 'L5.IT' ~ 'L3.5.', 
                                                                TRUE ~ new_cell_name)), 
               marker_df %>% mutate(new_cell_name = str_sub(subclass %>% make.names(), start = 1, end = 5)), 
               by = c('gene' ='gene', "new_cell_name" = "new_cell_name"))

# sort the data frame according to the following: finemap prob, snp id, cell type, and distance of gene from SNP
sorted_df = df %>% arrange(-finemap_posterior_probability, SNP, differentially_accessible_cell_type, distance)

# output sorted list of SNPs and cell types and genes < 100Mb away. filter for only genes that have enriched expression from 
# sonny's marker analysis or are the closest gene to the SNP (distance = 0)
sorted_df <- sorted_df %>% filter(distance < 100000, (avg_log2FC > 1 | distance == 0)) %>% 
  dplyr::select(SNP, accessible_region, finemap_posterior_probability, differentially_accessible_cell_type, accessibility_log_fold_change, gene, distance, avg_log2FC, pct.1, pct.2, expression, PsyOPS_score) 


#SNAREseq RNA differential expression
topSNPs <- read_csv('combined_notfilteredbypsyops.csv')
diff_exp_results <- readRDS("~/chromatin-accessibility/Adult/diff_exp_results.rds")
topSNPs_combined_RNAdiffexp <- left_join(topSNPs, diff_exp_results, by="gene")

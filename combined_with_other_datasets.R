library(data.table)
library(dplyr)
library(readr)
library(tidyverse)
library(janitor)

#Read in files
setwd("")
input <- read.table(file = 'PsyOPS_input.tsv', sep = '\t', header = TRUE)
output <- read.table(file = 'PsyOPS_output.tsv', sep = '\t', header = TRUE)
SZ <- read.csv(file="SZ_intersection_summary.csv")
gene_expression <- read_csv("new_CgG_avg_count.csv")

#Joins Psyops results with previously intersected SNPs
psyops_combined <- left_join(SZ, output, by=c("SNP"="lead_variant"))  %>%
  mutate(mapped_cell_types = case_when(!differentially_accessible_cell_type %in% c("L2.3.IT", "L5.IT", "L6.IT") ~ differentially_accessible_cell_type, 
                                       differentially_accessible_cell_type=="L2.3.IT" ~ "IT", 
                                       differentially_accessible_cell_type=="L5.IT" ~ "IT",
                                       differentially_accessible_cell_type=="L6.IT" ~ "IT"))%>%
  filter(P<5*10^-7 & finemap_posterior_probability>0.1 & PsyOPS_score>0.1)

#IT in gene expression sheet counts for 3 cell types in SNARE-seq dataset: L2.3.IT, L5.IT, L6.IT
#Joins gene expression with previously intersected SNPs
combined <- left_join(psyops_combined, gene_expression, by=c("gene"="X1"))%>%
  clean_names() %>%
  rename(Micro.PVM=microglia,
         Oligo=oligodendrocyte, 
         VIP=vip,
         OPC=opc, 
         Astro=astrocyte, 
         PVALB=pvalb, 
         LAMP5=lamp5,
         SST=sst, 
         IT=it) %>%
  pivot_longer(Oligo:l4_it,names_to = "differentially_expressed_cell_type", values_to = "expression") %>%
  filter(mapped_cell_types==differentially_expressed_cell_type & expression>1)


#Saves the list of genes in a line-separated text file
cat(unique(combined$gene), sep="\n", file="./gene_list.txt")

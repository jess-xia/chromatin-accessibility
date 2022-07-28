library(purrr)
library(plyr)


file_list <- list.files("/external/rprshnas01/kcni/jxia/chromatin-accessibility/MDD_parallel/parallel_results")
combined <- readRDS(paste("/external/rprshnas01/kcni/jxia/chromatin-accessibility/MDD_parallel/parallel_results/", file_list[1], sep=""))

for (i in 2:length(file_list)){
  combined <- Map(c, combined, readRDS(paste("/external/rprshnas01/kcni/jxia/chromatin-accessibility/MDD_parallel/parallel_results/", file_list[i], sep="")))
}

combined_df <- as.data.frame(do.call(cbind, combined)) %>% #Converts nested list to dataframe
  #Saves only rows that have a non-zero value
  filter(Astro!=0|Endo!=0|`L2-3 IT`!=0|`L5 ET`!=0| `L5 IT`!=0| `L5-6 NP`!=0| `L6 CT`!=0| `L6 IT`!=0| `L6 IT Car3`!=0| L6b!=0| LAMP5!=0| `Micro-PVM`!=0| Oligo!=0| OPC!=0| PVALB!=0| SNCG!=0| SST!=0| `SST CHODL`!=0| VIP!=0| VLMC!=0)
saveRDS(combined_df, file="MDD_intersected.rds")

library(purrr)
library(plyr)
library(tidyverse)

if(FALSE){ #Developmental, needs more work
file_list <- list.files("/external/rprshnas01/kcni/jxia/chromatin-accessibility/Developmental/SZ_parallel/parallel_results")
combined <- readRDS(paste("/external/rprshnas01/kcni/jxia/chromatin-accessibility/Developmental/SZ_parallel/parallel_results/", file_list[1], sep=""))

for (i in 2:length(file_list)){
  combined <- Map(c, combined, readRDS(paste("/external/rprshnas01/kcni/jxia/chromatin-accessibility/Developmental/SZ_parallel/parallel_results/", file_list[i], sep="")))
}

combined_df <- as.data.frame(do.call(cbind, combined)) %>% #Converts nested list to dataframe
  #Saves only rows that have a non-zero value
  filter(Astro!=0|Endo!=0|`L2-3 IT`!=0|`L5 ET`!=0| `L5 IT`!=0| `L5-6 NP`!=0| `L6 CT`!=0| `L6 IT`!=0| `L6 IT Car3`!=0| L6b!=0| LAMP5!=0| `Micro-PVM`!=0| Oligo!=0| OPC!=0| PVALB!=0| SNCG!=0| SST!=0| `SST CHODL`!=0| VIP!=0| VLMC!=0)
saveRDS(combined_df, file="SZ_intersected.rds")

write.csv(combined_df, "~/chromatin-accessibility/Developmental/SZ_intersected.csv")
}



# For Adult data
file_list <- list.files("/external/rprshnas01/kcni/jxia/chromatin-accessibility/Adult/SZ_parallel/parallel_results")
combined <- as.data.frame(readRDS(paste("/external/rprshnas01/kcni/jxia/chromatin-accessibility/Adult/SZ_parallel/parallel_results/", file_list[1], sep="")))

for (i in 2:length(file_list)){
  combined <- rbind(combined, as.data.frame(readRDS(paste("/external/rprshnas01/kcni/jxia/chromatin-accessibility/Adult/SZ_parallel/parallel_results/", file_list[i], sep=""))))
}

narrowed_down <- combined %>% 
  pivot_longer(Astro:VLMC,names_to = "cell_type", values_to = "accessibility") %>%
  filter(accessibility > 0)

saveRDS(narrowed_down, file="~/chromatin-accessibility/Adult/SZ_intersected_hg38.rds")
  
write.csv(narrowed_down, "~/Adult/SZ_intersected_hg38.csv")
  

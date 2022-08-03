library(readxl)
library(dplyr)

peak_range <- read_xlsx("/external/rprshnas01/kcni/jxia/chromatin-accessibility/Developmental/Developmental_DAR.xlsx", sheet=1) %>%
  mutate(chr=str_match(seqnames, "chr([A-Z0-9]{1,2})")[,2]) %>%
  subset(select=c("chr", "start", "end","peak_name"))
  

cell_type_specific_peaks <- read_xlsx("/external/rprshnas01/kcni/jxia/chromatin-accessibility/Developmental/Developmental_DAR.xlsx", sheet=3)
DAR_developmental <- left_join(cell_type_specific_peaks, peak_range, by="peak_name")

saveRDS(DAR_developmental, file="~/chromatin-accessibility/Developmental/Developmental_DAR.rds")

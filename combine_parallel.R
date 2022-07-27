library(purrr)

file_list <- list.files("/external/rprshnas01/kcni/jxia/chromatin-accessibility/SZ_parallel/parallel_results")
combined <- readRDS(paste("/external/rprshnas01/kcni/jxia/chromatin-accessibility/SZ_parallel/parallel_results/", file[1], sep=""))

for (i in 2:length(file_list)){
  combined <- Map(c, combined, readRDS(paste("/external/rprshnas01/kcni/jxia/chromatin-accessibility/SZ_parallel/parallel_results/", file_list[i], sep="")))
}
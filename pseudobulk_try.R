library(scran)
library(scuttle)
set.seed(10000)

sce <- mockSCE(ncells=1000)
sce$samples <- gl(8, 125) # Pretending we have 8 samples.

# Making up some clusters.
sce <- logNormCounts(sce)
clusters <- kmeans(t(logcounts(sce)), centers=3)$cluster

# Creating a set of pseudo-bulk profiles:
info <- DataFrame(sample=sce$samples, cluster=clusters)
pseudo <- sumCountsAcrossCells(sce, info)

# Making up an experimental design for our 8 samples.
pseudo$DRUG <- gl(2,4)[pseudo$sample]

# DGE analysis:
out <- pseudoBulkDGE(pseudo, 
                     label=pseudo$cluster,
                     condition=pseudo$DRUG,
                     design=~DRUG,
                     coef="DRUG2"
)
out[[1]]
metadata(out[[1]])$design


dge <- pseudoBulkDGE(x=RNA_counts, 
                     col.data=meta_data, 
                     label=meta_data$subclass)







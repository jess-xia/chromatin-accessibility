library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)
library(patchwork)
library(DESeq2)
library(textshape)

RNA_counts <- readRDS("/external/rprshnas01/netdata_kcni/stlab/ATACseq/Zhang_BICCN-H_20190523-20190611_huMOp_Final_RNA_Counts.RDS")
meta_data <- fread("/external/rprshnas01/netdata_kcni/stlab/ATACseq/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt")

#View the sparse matrix
#RNA_counts[1:50, 1:50]

## pseudobulk gene expression per cell-type
getPseudobulk <- function(mat, celltype) {
  mat.summary <- do.call(cbind, lapply(levels(celltype), function(ct) {
    cells <- names(celltype)[celltype==ct]
    pseudobulk <- rowSums(mat[, cells])
    return(pseudobulk)
  }))
  colnames(mat.summary) <- levels(celltype)
  return(mat.summary)
}

#Use sample_name column as rownames
meta_data <- column_to_rownames(meta_data, 'sample_name')

#Create object that maps samples to cell subtypes
celltype <- meta_data$subclass
levels(celltype) <- paste0(unique(meta_data$subclass))
names(celltype) <- rownames(meta_data)

#Checks if the sample names are identical & in the same order
#identical(colnames(RNA_counts), rownames(meta_data))

#Generate count matrix the sums up the counts for each gene for each cell type
counts <- getPseudobulk(RNA_counts, celltype)

#Sum up counts for each gene(row)
#counts_per_gene <- Matrix::rowSums(counts)



#Create Seurat object using sparse matrix and meta data
SNAREseq.data <- RNA_counts
SNAREseq <- CreateSeuratObject(counts = SNAREseq.data, project = "SNAREseq", min.cells = 3, min.features = 200, meta.data = meta_data)
Idents(object = SNAREseq) <- "subclass"
SNAREseq <- NormalizeData(SNAREseq, normalization.method = "LogNormalize", scale.factor = 10000)

#head(Idents(SNAREseq), 5)
#Idents(SNAREseq)

SNAREseq.markers <- FindAllMarkers(SNAREseq, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
diff_exp_results <- SNAREseq.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

saveRDS(diff_exp_results, file = "diff_exp_results.rds")


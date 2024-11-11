library(Seurat)
library(SeuratDisk)
library(SeuratData)

sciPlex_cds <- readRDS("GSM7056149_sciPlexGxE_2_preprocessed_cds.rds")

cds_A172 <- as.Seurat(sciPlex_cds$A172, data = NULL)
cds_T98G <- as.Seurat(sciPlex_cds$T98G, data = NULL)
cds_U87MG <- as.Seurat(sciPlex_cds$U87MG, data = NULL)

SaveH5Seurat(cds_A172, filename = "cds_A172.h5Seurat")
SaveH5Seurat(cds_T98G, filename = "cds_T98G.h5Seurat")
SaveH5Seurat(cds_U87MG, filename = "cds_U87MG.h5Seurat")

Convert("cds_A172.h5Seurat", dest = "h5ad")
Convert("cds_T98G.h5Seurat", dest = "h5ad")
Convert("cds_U87MG.h5Seurat", dest = "h5ad")

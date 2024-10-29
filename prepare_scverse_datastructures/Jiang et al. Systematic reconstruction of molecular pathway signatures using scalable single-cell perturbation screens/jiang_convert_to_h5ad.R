library(sceasy)
library(reticulate)
library(Seurat)

anndata <- reticulate::import('anndata')

ifnb <- readRDS(file = "Seurat_object_IFNB_Perturb_seq.rds")
ifng <- readRDS(file = "Seurat_object_IFNG_Perturb_seq.rds")
ins <- readRDS(file = "Seurat_object_INS_Perturb_seq.rds")
tgfb <- readRDS(file = "Seurat_object_TGFB_Perturb_seq.rds")
tnfa <- readRDS(file = "Seurat_object_TNFA_Perturb_seq.rds")

sceasy::convertFormat(ifnb, from="seurat", to="anndata", outFile='IFNB_Perturb_seq.h5ad')
sceasy::convertFormat(ifng, from="seurat", to="anndata", outFile='IFNG_Perturb_seq.h5ad')
sceasy::convertFormat(ins, from="seurat", to="anndata", outFile='INS_Perturb_seq.h5ad')
sceasy::convertFormat(tgfb, from="seurat", to="anndata", outFile='TGFB_Perturb_seq.h5ad')
sceasy::convertFormat(tnfa, from="seurat", to="anndata", outFile='TNFA_Perturb_seq.h5ad')

## Convert DIALOGUEs cell.type object into anndata

# Compile S4 objects into Seurat object
#devtools::install_github("livnatje/DIALOGUE")
library(DIALOGUE)
library(Seurat)

# Download sample dataset from: https://github.com/livnatje/DIALOGUE/tree/master/Data
ex <- readRDS('~/Downloads/test.example.rds')

# create seurat objects
a = Seurat::CreateSeuratObject(counts = ex$A@tpm, meta.data = ex$A@metadata)
b = Seurat::CreateSeuratObject(counts = ex$B@tpm, meta.data = ex$B@metadata)
c = Seurat::CreateSeuratObject(counts = ex$C@tpm, meta.data = ex$C@metadata)

# add subset identifier
a$subset = 'A'
b$subset = 'B'
c$subset = 'C'

# merge
example = merge(a, y=c(b, c))

# Translation from Seurat Object to Anndata using sceasy
# https://github.com/cellgeni/sceasy
library(reticulate)
library(sceasy)
library(SeuratData)

# Conda environment as described in sceasy description
## https://github.com/cellgeni/sceasy
use_condaenv(condaenv = 'sceasy', required = TRUE)

# conversion
filename <- paste('ex_DIALOGUE_sceasy', '.h5ad', sep = '')
sceasy::convertFormat(example, from="seurat", to="anndata", outFile=filename)

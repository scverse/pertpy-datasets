# Translationg from Seurat Object to Anndata using sceasy:
# https://github.com/cellgeni/sceasy
library(reticulate)
library(sceasy)
library(SeuratData)


# Download sample dataset using SeuratData.
#InstallData(ds = "thp1.eccite")


# Load object.
eccite <- LoadData(ds = "thp1.eccite")

# Conda environment as described in sceasy description
use_condaenv(condaenv = 'py36', required = TRUE)


# Conversion of assays into Anndata.
assays <- c('RNA', 'GDO', 'HTO', 'ADT')

for(a in assays){
  filename <- paste('eccite_sceasy_', a, '.h5ad', sep = '')
  sceasy::convertFormat(eccite, from="seurat", to="anndata",
                        outFile=filename, assay=a)
}

The link to download the dataset can be found here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7056149

The downloaded directory contains a single RDS file that has to be loaded into your R environment. The RDS file contains a list of Seurat objects, each corresponding to a different cell line.

Then the Seurat objects are converted to a h5ad files and concatenated into a single h5ad file, which will be around 35GB if not compressed.

#!/usr/bin/env Rscript

require(Seurat)
require(SeuratDisk)

source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

scrnaID=arguments[1]

scrna_data=readRDS(paste0("scrna/",scrnaID,".rds"))


DefaultAssay(scrna_data) <- "RNA"

scrna_data@assays$SCT <- NULL


UpdateSeuratObject(scrna_data) -> scrna_data

SaveH5Seurat(scrna_data,paste0(scrnaID,".h5Seurat"))
SeuratDisk::Convert(paste0(scrnaID,".h5Seurat"), dest = "h5ad")






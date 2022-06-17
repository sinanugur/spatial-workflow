#!/usr/bin/env Rscript

require(Seurat)
require(SeuratDisk)

source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

scrnaID=arguments[1]

scrna_data=readRDS(paste0("scrna/",scrnaID,".rds"))

UpdateSeuratObject(scrna_data) -> scrna_data
DefaultAssay(scrna_data) <- "RNA"
scrna_data[["SCT"]] <- NULL

DietSeurat(scrna_data) -> scrna_data


SaveH5Seurat(scrna_data,paste0("scrna/",scrnaID,".h5Seurat"))
SeuratDisk::Convert(paste0("scrna/",scrnaID,".h5Seurat"), dest = "h5ad")







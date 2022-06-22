#!/usr/bin/env Rscript

require(Seurat)
require(tidyverse)
require(viridis)



source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

sampleID=arguments[1]
scrnaID=arguments[2]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))
dwls_data=readRDS(paste0("DWLS_assay/",scrnaID,"/",sampleID,".rds"))


function_image_fixer(Spatial_Data) -> Spatial_Data



Spatial_Data[["DWLS"]] <- dwls_data


DefaultAssay(Spatial_Data) <- "DWLS"

cell_types_all=rownames(Spatial_Data)

wp=seurat_plotting()

ggsave(paste0("results/",sampleID,"/deconvolution/dwls/",sampleID,"-",scrnaID,"-dwls.pdf"),wp,height=18,width=8)

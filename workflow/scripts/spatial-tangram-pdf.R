#!/usr/bin/env Rscript

require(Seurat)
require(tidyverse)
require(viridis)


source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

sampleID=arguments[1]
scrnaID=arguments[2]


Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))
tangram_csv=read.csv(paste0("tangram/",scrnaID,"/",sampleID,".csv"))



function_image_fixer(Spatial_Data) -> Spatial_Data

Spatial_Data[["predictions"]] <- CreateAssayObject(tangram_csv  %>% column_to_rownames("X") %>% t())


DefaultAssay(Spatial_Data) <- "predictions"


cell_types_all=tangram_csv  %>% column_to_rownames("X") %>% colnames()


wp=seurat_plotting()


ggsave(paste0("results/",sampleID,"/deconvolution/tangram/",sampleID,"-",scrnaID,"-tangram.pdf"),wp,height=18,width=8)

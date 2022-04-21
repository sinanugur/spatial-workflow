#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)
source("workflow/scripts/spatial-functions.R")

arguments=commandArgs(TRUE)

sampleID=arguments[1]
res=arguments[2]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))

function_image_fixer(Spatial_Data) -> Spatial_Data

p1 <- DimPlot(Spatial_Data, reduction = "umap", label = TRUE,label.size = 10,group.by = paste0("SCT_snn_res.",res))
p2 <- SpatialDimPlot(Spatial_Data, label = TRUE, label.size = 6,group.by = paste0("SCT_snn_res.",res),images=paste0("image"))

ggsave(plot =p1+p2,filename=paste0(sampleID,"/resolution-",res,"/",sampleID,".umap.spatial",".pdf"),width=13,height=7)



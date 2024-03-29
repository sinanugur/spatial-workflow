#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)
require(patchwork)

source("workflow/scripts/spatial-functions.R")

arguments=commandArgs(TRUE)

sampleID=arguments[1]
res=arguments[2]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))

function_image_fixer(Spatial_Data) -> Spatial_Data

p1 <- DimPlot(Spatial_Data, reduction = "umap", label = TRUE,label.size = 10,group.by = paste0("SCT_snn_res.",res)) 
p2 <- SpatialDimPlot(Spatial_Data, label = TRUE, label.size = 6,group.by = paste0("SCT_snn_res.",res),images=paste0("image")) 

suppressWarnings(wrap_plots(p1,p2,ncol=2) -> wp)


ggsave(plot =wp,filename=paste0("results/",sampleID,"/resolution-",res,"/",sampleID,".umap.spatial",".pdf"),width=13,height=7)



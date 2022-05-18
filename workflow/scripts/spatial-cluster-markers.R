#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)


arguments=commandArgs(TRUE)

sampleID=arguments[1]
res=arguments[2]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))

SCT_=paste0("SCT_snn_res.",res)


Idents(object = Spatial_Data) <- Spatial_Data@meta.data[[SCT_]]


all_markers=FindAllMarkers(Spatial_Data)

openxlsx::write.xlsx(all_markers,file=paste0("results/",sampleID,"/resolution-",res,"/",sampleID,".all-markers-forAllClusters",".xlsx"))

openxlsx::write.xlsx(all_markers %>% filter(avg_log2FC > 0),file=paste0("results/",sampleID,"/resolution-",res,"/",sampleID,".positive-markers-forAllClusters",".xlsx"))




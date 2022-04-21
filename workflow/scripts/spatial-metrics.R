#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)


arguments=commandArgs(TRUE)

sampleID=arguments[1]
res=arguments[2]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))





SCT_=paste0("SCT_snn_res.",res)

metrics=table(Spatial_Data@meta.data[[SCT_]], Spatial_Data@meta.data$orig.ident)

openxlsx::write.xlsx(metrics %>% as.data.frame() %>% select(sampleID=1),file=paste0(sampleID,"/resolution-",res,"/",sampleID,".number-of-cells-per-cluster",".xlsx"))

#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)


arguments=commandArgs(TRUE)

sampleID=arguments[1]


Spatial_Data=readRDS(paste0(sampleID,".rds"))

Spatial_Data <- FindClusters(Spatial_Data, verbose = FALSE,resolution = seq(0.1,2.5,0.1))


saveRDS(Spatial_Data,file = paste0("../",sampleID,".rds"))
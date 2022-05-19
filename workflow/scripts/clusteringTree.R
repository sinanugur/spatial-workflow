#!/usr/bin/env Rscript

library(Seurat)
require(tidyverse)
require(clustree)
arguments=commandArgs(TRUE)

sampleID=arguments[1]


Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))

Spatial_Data <- FindClusters(Spatial_Data, verbose = FALSE,resolution = seq(0.1,2.5,0.1))

clustree(Spatial_Data) -> p1

ggsave(paste0("results/",sampleID,"/clusteringTree/clusteringTree-",sampleID,".pdf"),p1,width=8,height=15)
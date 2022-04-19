#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)


arguments=commandArgs(TRUE)

sampleID=arguments[1]

Spatial_Df=readRDS(paste0("rds/",sampleID,".rds"))


vp=VlnPlot(Spatial_Df, features = c( "nCount_Spatial","nFeature_Spatial"),group.by = "orig.ident", pt.size = 0.1) + NoLegend() 

ggsave(vp,paste0("pdf/",sampleID,".violin",".pdf"))


p1 <- DimPlot(Spatial_Df, reduction = "umap", label = TRUE,label.size = 10)
p2 <- SpatialDimPlot(Spatial_Df, label = TRUE, label.size = 6)

ggsave(p1+p2,paste0("pdf/",sampleID,".umap",".pdf"))
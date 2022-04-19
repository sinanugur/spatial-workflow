#!/usr/bin/env Rscript

library(Seurat)
require(tidyverse)

source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

sampleID=arguments[1]
slice=arguments[2]
datadir=c(paste0(sampleID,"/outs/"))

Spatial_Df=function_analyse_spatial(datadir = datadir,slice = slice)
saveRDS(Spatial_Df,file = paste0("rds/",sampleID,".rds"))


#vp=VlnPlot(Spatial_Df, features = c( "nCount_Spatial","nFeature_Spatial"),group.by = "orig.ident", pt.size = 0.1) + NoLegend() 

#ggsave(vp,paste0("pdf/",sampleID,".violin",".pdf"))


#p1 <- DimPlot(tuvas_131, reduction = "umap", label = TRUE,label.size = 10)
#p2 <- SpatialDimPlot(tuvas_131, label = TRUE, label.size = 6)

#ggsave(p1+p2,paste0("pdf/",sampleID,".umap",".pdf"))
#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)


arguments=commandArgs(TRUE)

sampleID=arguments[1]
res=arguments[2]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))

IMAGE=Read10X_Image(image.dir=paste0("data/",sampleID,"/outs/spatial"),image.name="grayscale.png")

Spatial_Data@images$"image" <- IMAGE

Spatial_Data@images$"image"@assay <- "Spatial"
Spatial_Data@images$"image"@assay <- "Spatial"

Spatial_Data@images$"image"@key <- paste0("image","_")

p1 <- DimPlot(Spatial_Data, reduction = "umap", label = TRUE,label.size = 10,group.by = paste0("SCT_snn_res.",res))
p2 <- SpatialDimPlot(Spatial_Data, label = TRUE, label.size = 6,group.by = paste0("SCT_snn_res.",res),images=paste0("image"))

ggsave(plot =p1+p2,filename=paste0(sampleID,"/resolution-",res,"/",sampleID,".umap.spatial",".pdf"),width=13,height=7)



#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)
require(patchwork)


arguments=commandArgs(TRUE)

sampleID=arguments[1]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))



IMAGE=Read10X_Image(image.dir=paste0("data/",sampleID,"/outs/spatial"),image.name="grayscale.png")



Spatial_Data@images$"image" <- IMAGE

Spatial_Data@images$"image"@assay <- "Spatial"
Spatial_Data@images$"image"@assay <- "Spatial"

Spatial_Data@images$"image"@key <- paste0("image","_")


plot1 <- VlnPlot(Spatial_Data, features = "nCount_Spatial", pt.size = 0.1,assay="Spatial",group.by="orig.ident") + NoLegend()
plot2 <- SpatialFeaturePlot(Spatial_Data, features = "nCount_Spatial",images="image") + theme(legend.position = "right")

wrap_plots(plot1, plot2,ncol=2)

ggsave(filename=paste0(sampleID,"/technicals/",sampleID,".n_counts.pdf"),width=13,height=7)
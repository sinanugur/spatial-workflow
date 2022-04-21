#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)


arguments=commandArgs(TRUE)

sampleID=arguments[1]


Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))

IMAGE=Read10X_Image(image.dir=paste0("data/",sampleID,"/outs/spatial"),image.name="grayscale.png")

Spatial_Data@images$"image" <- IMAGE

Spatial_Data@images$"image"@assay <- "Spatial"
Spatial_Data@images$"image"@assay <- "Spatial"

Spatial_Data@images$"image"@key <- paste0("image","_")


markers=SpatiallyVariableFeatures(Spatial_Data, 
selection.method = "markvariogram") 

openxlsx::write.xlsx(markers %>% as.data.frame() %>% select(gene=1),file=paste0(sampleID,"/",sampleID,".spatial_markers",".xlsx"))

for (i in markers) {

SpatialFeaturePlot(Spatial_Data, features = i, ncol = 1, alpha = c(0.1, 1),images=paste0("image"))
ggsave(paste0(sampleID,"/","spatial-markers/",i,".pdf"))

}
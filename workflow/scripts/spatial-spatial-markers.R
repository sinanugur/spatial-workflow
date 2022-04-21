#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)
source("workflow/scripts/spatial-functions.R")

arguments=commandArgs(TRUE)

sampleID=arguments[1]


Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))

function_image_fixer(Spatial_Data) -> Spatial_Data


markers=SpatiallyVariableFeatures(Spatial_Data, 
selection.method = "markvariogram") 

openxlsx::write.xlsx(markers %>% as.data.frame() %>% select(gene=1),file=paste0(sampleID,"/spatial-markers/",sampleID,".spatial_markers",".xlsx"))

for (i in markers) {

SpatialFeaturePlot(Spatial_Data, features = i, ncol = 1, alpha = c(0.1, 1),images=paste0("image"))
ggsave(paste0(sampleID,"/","spatial-markers/plots/",i,".pdf"))

}
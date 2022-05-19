#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)
source("workflow/scripts/spatial-functions.R")

arguments=commandArgs(TRUE)

sampleID=arguments[1]


Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))

function_image_fixer(Spatial_Data) -> Spatial_Data



domanska_muscularis=data.frame(gene=c("FCER1A","CDC1C","CLEC10A","CCL3L1","CCL3","CCL4L2","MT1X","MT1E","CTSL","RGS1","FOS","APOE","DNASE1L3","MMP9","LYZ","AREG","EREG","CCL20","S100A9","S100A8","EREG","FCN1","VCAN","LYZ","HSPA1A","HSPA6","HSPA1B","LYVE1","MARCO","COLEC12"),group=c(3,3,3,5,5,5,6,6,6,7,7,7,4,4,4,1,1,1,0,0,0,2,2,2,9,9,9,11,11,11))

domanska_mucosas=data.frame(gene=c("LYVE1","F13A1","FOLR2","SELENOP","APOE","SLC40A1","C1QB","DAB2","PDK4","SPP1","ACP5","CD9","FCER1A","CD1C","CLEC10A","HSPA6","DNAJB1","HSPA1B","S100A8","S100A9","S100A12","EREG","G0S2","FCN1","CCL20","IL1B","IL23A","CXCL10","CXCL9","GBP1"),group=c(12,12,12,11,11,11,10,10,10,9,9,9,8,8,8,5,5,5,0,0,0,3,3,3,2,2,2,6,6,6))

domanska_markers=bind_rows(domanska_mucosas,domanska_muscularis) %>% distinct(gene) %>% pull()



for (i in domanska_markers) {

try({SpatialFeaturePlot(Spatial_Data, features = i, ncol = 1, alpha = c(0.1, 1),images=paste0("image"),pt.size.factor=1.1) + scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlGn")))
ggsave(paste0("results/",sampleID,"/","selected-markers/plots/",i,".pdf"))})

}
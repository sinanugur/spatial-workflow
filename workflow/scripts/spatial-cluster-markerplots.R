#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)
require(patchwork)
source("workflow/scripts/spatial-functions.R")

arguments=commandArgs(TRUE)

sampleID=arguments[1]
res=arguments[2]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))

function_image_fixer(Spatial_Data) -> Spatial_Data

SCT_=paste0("SCT_snn_res.",res)


Idents(object = Spatial_Data) <- Spatial_Data@meta.data[[SCT_]]


Positive_Features=openxlsx::read.xlsx(paste0(sampleID,"/resolution-",res,"/",sampleID,".positive-markers-forAllClusters.xlsx")) %>% select(cluster,gene) 


for(d in (Positive_Features %>% distinct(cluster) %>% pull())) {

    dir.create(paste0(sampleID,"/resolution-",res,"/markers/","cluster",d,"/"),recursive=TRUE)
}

options(warn=-1)
suppressMessages(for (i in 1:nrow(Positive_Features)) {

    gene=Positive_Features[i,]$gene
    cluster=Positive_Features[i,]$cluster

p1 <- FeaturePlot(Spatial_Data, reduction = "umap", features=gene) + scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlGn")))
p2 <- SpatialFeaturePlot(Spatial_Data, features=gene,images=paste0("image"),alpha=c(0.1,1),pt.size.factor=1.2) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlGn")))
p3 <- DotPlot(Spatial_Data, features=gene)
p4 <- VlnPlot(Spatial_Data,features=gene)

suppressWarnings(wrap_plots(p1,p2,p3,p4,ncol=2) -> wp)

ggsave(paste0(sampleID,"/resolution-",res,"/markers/","cluster",cluster,"/",gene,".pdf"),wp,height=14,width=14)

})






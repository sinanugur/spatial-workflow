#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)
require(patchwork)


arguments=commandArgs(TRUE)

sampleID=arguments[1]
res=arguments[2]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))

IMAGE=Read10X_Image(image.dir=paste0("data/",sampleID,"/outs/spatial"),image.name="grayscale.png")



Spatial_Data@images$"image" <- IMAGE

Spatial_Data@images$"image"@assay <- "Spatial"
Spatial_Data@images$"image"@assay <- "Spatial"

Spatial_Data@images$"image"@key <- paste0("image","_")

SCT_=paste0("SCT_snn_res.",res)


Idents(object = Spatial_Data) <- Spatial_Data@meta.data[[SCT_]]


Positive_Features=openxlsx::read.xlsx(paste0(sampleID,"/resolution-",res,"/",sampleID,".positive-markers-forAllClusters.xlsx")) %>% select(cluster,gene)


for(d in (Positive_Features %>% distinct(cluster))) {

    dir.create(paste0(sampleID,"/resolution-",res,"/markers/","cluster",d,"/"))
}

for (i in 1:nrow(Positive_Features)) {

    gene=Positive_Features[i,]$gene
    cluster=Positive_Features[i,]$cluster

p1 <- FeaturePlot(Spatial_Data, reduction = "umap", features=gene) + scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
p2 <- SpatialFeaturePlot(Spatial_Data, features=gene,images=paste0("image"),alpha=c(0.8,1),pt.size.factor=1.2) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
p3 <- DotPlot(Spatial_Data, features=gene)
p4 <- VlnPlot(Spatial_Data,features=gene)

wrap_plots(p1,p2,p3,p4,ncol=2)

ggsave(paste0(sampleID,"/resolution-",res,"/markers/","cluster",cluster,"/",gene,".pdf"),height=14,width=14)

}






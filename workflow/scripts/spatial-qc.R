#!/usr/bin/env Rscript
library(Seurat)
require(tidyverse)
require(patchwork)
source("workflow/scripts/spatial-functions.R")

arguments=commandArgs(TRUE)

sampleID=arguments[1]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))

Spatial_Data[["percent.mt"]] <- PercentageFeatureSet(Spatial_Data, pattern = "^MT-")
Spatial_Data[["percent.rp"]] <- PercentageFeatureSet(Spatial_Data, pattern = "^RP[SL]")



function_image_fixer(Spatial_Data) -> Spatial_Data


plot1 <- VlnPlot(Spatial_Data, features = "nCount_Spatial", pt.size = 0.1,assay="Spatial",group.by="orig.ident") + NoLegend()
plot2 <- SpatialFeaturePlot(Spatial_Data, features = "nCount_Spatial",images="image") + theme(legend.position = "right")

wrap_plots(plot1, plot2,ncol=2)

ggsave(filename=paste0("results/",sampleID,"/technicals/",sampleID,".n_counts.pdf"),width=13,height=7)


Spatial_Data <- NormalizeData(Spatial_Data, verbose = FALSE, assay = "Spatial")

Spatial_Data <- GroupCorrelation(Spatial_Data, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
Spatial_Data <- GroupCorrelation(Spatial_Data, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)

p1 <- GroupCorrelationPlot(Spatial_Data, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") + theme(plot.title = element_text(hjust = 0.5))
p2 <- GroupCorrelationPlot(Spatial_Data, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") + theme(plot.title = element_text(hjust = 0.5))


wrap_plots(p1, p2,ncol=2)

ggsave(filename=paste0("results/",sampleID,"/technicals/",sampleID,".normalization.pdf"),width=13,height=7)
#!/usr/bin/env Rscript

require(Seurat)
require(tidyverse)
require(viridis)


source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

sampleID=arguments[1]
scrnaID=arguments[2]


Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))
tangram_csv=read.csv(paste0("tangram/",scrnaID,"/",sampleID,".csv"))





Spatial_Data[["predictions"]] <- CreateAssayObject(tangram_csv  %>% column_to_rownames("X") %>% t())


DefaultAssay(Spatial_Data) <- "predictions"


cell_types_all=tangram_csv  %>% column_to_rownames("X") %>% colnames()


wp=Seurat::SpatialFeaturePlot(
  object = Spatial_Data,
  features = cell_types_all,alpha = c(0.7, 1),pt.size.factor = 1.5,ncol=2,images=paste0("image")) & 
  scale_fill_continuous(type = "viridis",labels=mylabels) &   
  theme(legend.title = element_text(size=4),legend.key.size = unit(0.2,"cm"),
  legend.text = element_text(size=3),legend.margin=margin(t = 0,b = 0, unit='cm'),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))


ggsave(paste0("results/",sampleID,"/deconvolution/tangram/",sampleID,"-",scrnaID,"-tangram.pdf"),wp,height=18,width=8)

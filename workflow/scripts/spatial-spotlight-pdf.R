#!/usr/bin/env Rscript

library(Seurat)
require(tidyverse)
library(SPOTlight)


source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

sampleID=arguments[1]
scrnaID=arguments[2]

Spatial_Data=readRDS(paste0("rds_decon/",scrnaID,"/",sampleID,".rds"))
scrna_data=readRDS(paste0("scrna/",scrnaID,".rds"))


cell_types_all=Idents(scrna_data) %>% unique() %>% as.character() %>% make.names()

function_image_fixer(Spatial_Data) -> Spatial_Data


wp=Seurat::SpatialFeaturePlot(
  object = Spatial_Data,
  features = cell_types_all,
  alpha = c(0.1, 1),pt.size.factor = 1,ncol=2,images=paste0("image")) & scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlGn"))) & theme(legend.title = element_text(size=4),legend.key.size = unit(0.2,"cm"),legend.text = element_text(size=3),legend.margin=margin(t = 0,b = 0, unit='cm'),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))


ggsave(paste0("results/",sampleID,"/deconvolution/spotlight/",sampleID,"-",scrnaID,"-spotlight.pdf"),wp,height=18,width=6)

#!/usr/bin/env Rscript

library(Seurat)
require(tidyverse)
library(SPOTlight)
require(viridis)


source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

sampleID=arguments[1]
scrnaID=arguments[2]

Spatial_Data=readRDS(paste0("rds_decon/",scrnaID,"/",sampleID,".rds"))
scrna_data=readRDS(paste0("scrna/",scrnaID,".rds"))


cell_types_all=Idents(scrna_data) %>% unique() %>% as.character() %>% make.names()

function_image_fixer(Spatial_Data) -> Spatial_Data


wp=seurat_plotting()

ggsave(paste0("results/",sampleID,"/deconvolution/spotlight/",sampleID,"-",scrnaID,"-spotlight.pdf"),wp,height=18,width=6)

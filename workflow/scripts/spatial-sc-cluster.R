#!/usr/bin/env Rscript

library(Seurat)
require(tidyverse)
library(SPOTlight)


source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

scrnaID=arguments[1]

scrna_data=readRDS(paste0("scrna/",scrnaID,".rds"))





cluster_markers_all <- Seurat::FindAllMarkers(object = scrna_data, 
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)


openxlsx::write.xlsx(cluster_markers_all,file=paste0("scrna/",scrnaID,".cluster_markers.xlsx"))

#!/usr/bin/env Rscript

library(Seurat)
require(tidyverse)

source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

sampleID=arguments[1]
slice=make.names(arguments[1])
datadir=c(paste0("data/",sampleID,"/outs/"))

Spatial_Data=function_analyse_spatial(datadir = datadir,slice = slice)
saveRDS(Spatial_Data,file = paste0("rds/",sampleID,".rds"))



#!/usr/bin/env Rscript

library(Seurat)
require(tidyverse)
library(SPOTlight)


source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

sampleID=arguments[1]
scrnaID=arguments[2]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))
scrna_data=readRDS(paste0("scrna/",scrnaID,".rds"))


function_image_fixer(Spatial_Data) -> Spatial_Data



openxlsx::read.xlsx(paste0("scrna/",scrnaID,".cluster_markers.xlsx")) -> cluster_markers_all

function_spotlight=function(Seurat_Object,scrna_rds) {
set.seed(123)

spotlight_ls <- spotlight_deconvolution(
  se_sc = scrna_rds,
  counts_spatial = Seurat_Object@assays$Spatial@counts,
  clust_vr = "seurat_clusters", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 3000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
  )
nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
rownames(decon_mtrx) <- colnames(Seurat_Object)

decon_df <- decon_mtrx %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

Seurat_Object@meta.data <- Seurat_Object@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

return(Seurat_Object)
}

function_spotlight(Spatial_Data,scrna_data) -> Spatial_Data

saveRDS(Spatial_Data,file = paste0("rds_decon/",scrnaID,"/",sampleID,".rds"))




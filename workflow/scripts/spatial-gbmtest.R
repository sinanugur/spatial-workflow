#!/usr/bin/env Rscript

require(Seurat)
require(tidyverse)
#require(caret)


source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

sampleID=arguments[1]
modelID=arguments[2]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))
model=readRDS(paste0("models/",modelID,".rds"))





function_add_missing=function(df,vector) {
  
  for (i in vector) {
    df <- df %>% dplyr::mutate("{i}":=0)
    
  }
  return(df)
}

function_image_fixer(Spatial_Data) -> Spatial_Data


function_gbm_test=function(Seurat_Object) {
  GetAssayData(Seurat_Object,slot = "scale.data") %>% as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("cell") %>% janitor::clean_names()-> testdf2

  
  testdf2 <-function_add_missing(testdf2,setdiff(c(model$coefnames,"cell"),colnames(testdf2)))
  
  predict(model,testdf2,type="prob") %>% janitor::clean_names() %>% bind_cols(testdf2 %>% select(barcodes=cell)) -> gbm_decon_df


cell_types_all=gbm_decon_df %>% select(-barcodes) %>% colnames()

Seurat_Object@meta.data <- Seurat_Object@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(gbm_decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")
  
wp=Seurat::SpatialFeaturePlot(
  object = Seurat_Object,
  features = cell_types_all,
  alpha = c(0.1, 1),pt.size.factor = 1,ncol=2,images=paste0("image")) & scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlGn")),limits=c(0,1)) 


ggsave(paste0(sampleID,"/deconvolution/gbm/",sampleID,"-",modelID,"-gbm.pdf"),wp,height=20,width=7)

return(Seurat_Object)
}

function_gbm_test(Spatial_Data) -> Spatial_Data





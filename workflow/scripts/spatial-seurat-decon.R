#!/usr/bin/env Rscript

require(Seurat)
require(tidyverse)
require(viridis)



params=list(k.anchor=23,k.score=5,k.filter=100,n.trees=100)


source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

sampleID=arguments[1]
scrnaID=arguments[2]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))
scrna_data=readRDS(paste0("scrna/",scrnaID,".rds"))


function_image_fixer(Spatial_Data) -> Spatial_Data


function_decon_seurat = function(reference,query,anc){

anchors <- FindTransferAnchors(reference = reference, query = query, normalization.method = "SCT",
k.anchor = anc,
k.score = params$k.score,
k.filter=params$k.filter,
n.trees = params$n.trees)

predictions.assay <- TransferData(anchorset = anchors, refdata = reference$seurat_clusters, prediction.assay = TRUE)
query[["predictions"]] <- predictions.assay


return(query)
}



for (i in c(params$k.anchor,20,15,10,5)) {

try({
  function_decon_seurat(reference=scrna_data,query=Spatial_Data,anc=i) -> Spatial_Data
  break
  }
  )


}

DefaultAssay(Spatial_Data) <- "predictions"

cell_types_all=Idents(scrna_data) %>% unique() %>% as.character()

wp=seurat_plotting()

ggsave(paste0("results/",sampleID,"/deconvolution/seurat/",sampleID,"-",scrnaID,"-seurat.pdf"),wp,height=18,width=8)

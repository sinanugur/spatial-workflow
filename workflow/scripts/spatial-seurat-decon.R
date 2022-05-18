#!/usr/bin/env Rscript

library(Seurat)
require(tidyverse)


source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

sampleID=arguments[1]
scrnaID=arguments[2]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))
scrna_data=readRDS(paste0("scrna/",scrnaID,".rds"))


function_image_fixer(Spatial_Data) -> Spatial_Data


function_decon_seurat = function(reference,query,k.weight=30){

anchors <- FindTransferAnchors(reference = reference, query = query, normalization.method = "SCT")

predictions.assay <- TransferData(anchorset = anchors, refdata = reference$seurat_clusters, prediction.assay = TRUE,
    weight.reduction = query[["pca"]], dims = 1:30,k.weight=k.weight)
query[["predictions"]] <- predictions.assay


return(query)
}




for (i in c(30,20)) {

try({
  function_decon_seurat(reference=scrna_data,query=Spatial_Data,k.weight=i) -> Spatial_Data
  break
  }
  )


}


DefaultAssay(Spatial_Data) <- "predictions"

cell_types_all=Idents(scrna_data) %>% unique() %>% as.character()

wp=Seurat::SpatialFeaturePlot(
  object = Spatial_Data,
  features = cell_types_all,
  alpha = c(0.1, 1),pt.size.factor = 1,ncol=2,images=paste0("image")) & scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlGn"))) & theme(legend.title = element_text(size=4),legend.key.size = unit(0.2,"cm"),legend.text = element_text(size=3),legend.margin=margin(t = 0,b = 0, unit='cm'),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))


ggsave(paste0("results/",sampleID,"/deconvolution/seurat/",sampleID,"-",scrnaID,"-seurat.pdf"),wp,height=18,width=8,useDingbats=TRUE)
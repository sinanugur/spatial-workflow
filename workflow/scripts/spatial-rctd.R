#!/usr/bin/env Rscript

require(Seurat)
require(tidyverse)
require(spacexr)


source("workflow/scripts/spatial-functions.R")
arguments=commandArgs(TRUE)

sampleID=arguments[1]
scrnaID=arguments[2]

Spatial_Data=readRDS(paste0("rds/",sampleID,".rds"))
scrna_data=readRDS(paste0("scrna/",scrnaID,".rds"))


cell_types_all=Idents(scrna_data) %>% unique() %>% as.character() %>% make.names()

function_image_fixer(Spatial_Data) -> Spatial_Data


counts <- data.frame(scrna_data@assays$RNA@counts)
colnames(counts) <- colnames(scrna_data)
meta_data <- data.frame(scrna_data@meta.data)
cell_types <- meta_data$seurat_clusters %>% make.names()
names(cell_types) <- rownames(scrna_data@meta.data)
cell_types <- as.factor(cell_types)
nUMI_df <- data.frame(colSums(scrna_data@assays$RNA@counts))
nUMI <- nUMI_df$colSums.scrna_data.assays.RNA
names(nUMI) <- rownames(nUMI_df)



reference <- Reference(counts, cell_types, nUMI)




coords <- data.frame(colnames(Spatial_Data))
colnames(coords) <- 'barcodes'
coords$xcoord <- seq_along(colnames(Spatial_Data))
coords$ycoord <- seq_along(colnames(Spatial_Data))
counts <- data.frame(Spatial_Data@assays$Spatial@counts) # load in counts matrix
colnames(counts) <- colnames(Spatial_Data)
coords <- data.frame(colnames(Spatial_Data))
colnames(coords) <- 'barcodes'
coords$xcoord <- seq_along(colnames(Spatial_Data))
coords$ycoord <- seq_along(colnames(Spatial_Data))
rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI

puck<- SpatialRNA(coords = coords,counts = counts,nUMI = nUMI)

myRCTD <- create.RCTD(puck,reference,max_cores = 5,UMI_min = 20)
myRCTD <- run.RCTD(myRCTD,doublet_mode = "full")


Spatial_Data@meta.data <- Spatial_Data@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(myRCTD@results$weights %>% as.data.frame() %>% rownames_to_column("barcodes"), by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")


saveRDS(Spatial_Data,file = paste0("rds_rctd/",scrnaID,"/",sampleID,".rds"))
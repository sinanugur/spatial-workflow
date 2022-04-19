
function_pca_dimensions=function(Spatial_Data){
  
  pct <- Stdev(object = Spatial_Data, reduction = "pca") / sum(Stdev(object = Spatial_Data, reduction = "pca")) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cum > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),  decreasing = T)[1] + 1 
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
dimensionReduction <- min(co1, co2) 
  
}

function_read_spatial=function(datadir,filename="filtered_feature_bc_matrix.h5",assay="Spatial",slice="slice1"){
  
  Spatial_Data <- Load10X_Spatial(data.dir = datadir, 
                                filename = filename,
                                assay = "Spatial",
                                slice = slice)

  return(Spatial_Data)
}

function_spatial_features=function(Spatial_Data) {
  
  Spatial_Data <- FindSpatiallyVariableFeatures(Spatial_Data, assay = "SCT", features = VariableFeatures(Spatial_Data)[1:1000],selection.method = "markvariogram")
  
}
  
function_analyse_spatial=function(datadir,filename="filtered_feature_bc_matrix.h5",assay="Spatial",slice="slice1",calculatespatialfeatures=T){
  

Spatial_Data=function_read_spatial(datadir=datadir,filename=filename,assay=assay,slice=slice)

  #transformation
Spatial_Data <- SCTransform(Spatial_Data, assay = "Spatial", verbose = FALSE)
  
  


Spatial_Data <- RunPCA(Spatial_Data, assay = "SCT", verbose = FALSE)
dimensionReduction=function_pca_dimensions(Spatial_Data)

Spatial_Data <- FindNeighbors(Spatial_Data, reduction = "pca", dims = 1:dimensionReduction)
Spatial_Data <- FindClusters(Spatial_Data, verbose = FALSE)
Spatial_Data <- RunUMAP(Spatial_Data, reduction = "pca", dims = 1:dimensionReduction)

  if(isTRUE(calculatespatialfeatures== TRUE)) {
      Spatial_Data <- function_spatial_features(Spatial_Data)

      }
}


##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(300000)

##### Load Packages #####
  #### Basic installation ####
  Package.set <- c("tidyverse","Seurat","ggplot2","ggpmisc",
                   "stringr","magrittr","dplyr")
  ## Check whether the installation of those packages is required from basic
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)


##### Function setting #####
  ## Call function
  source("FUN_Beautify_ggplot.R")

#### Load data #####
  load("SeuratObject_CDS_PRJCA001063.RData")


##### CombineSeuObj #####
  #### Re-dimension reduction ####
  # scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj, selection.method = "vst", nfeatures = 2000)
  scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj)
  # # Run the standard workflow for visualization and clustering
  # scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
  scRNA.SeuObj <- RunPCA(scRNA.SeuObj, npcs = 200, verbose = FALSE)
  # scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:160,n.neighbors = 20,min.dist = 0.3)
  scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:200)
  scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)

  ##### Plot #####
  scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:160,n.neighbors = 20, min.dist = 0.3)
  FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))
  DimPlot(scRNA.SeuObj, reduction = "umap")
  DimPlot(scRNA.SeuObj, reduction = "umap",label = T)
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type")

##### Save RData #####
  save.image("SeuratObject_CDS_PRJCA001063.RData")




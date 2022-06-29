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


##### Cell type ReAnnotation*  #####

  DimPlot(scRNA.SeuObj, label = TRUE)

  Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["leiden"]]
  DimPlot(scRNA.SeuObj, label = TRUE)

  Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["seurat_clusters"]]
  DimPlot(scRNA.SeuObj, label = TRUE)

  scRNA.SeuObj <- RenameIdents(scRNA.SeuObj, `0` = "CD4+T", `1` = "B", `2` = "Mac3",
                               `3` = "Neu", `4` = "CD8+T", `5` = "CD8+T", `6` = "Mac2", `7` = "CD4+T",
                               `8` = "NK", `9` = "Neu",`10` = "Mast1", `11` = "T", `12` = "Ery", `13` = "Mac1",
                               `14` = "B", `15` = "B", `16` = "Mast2", `17` = "Mac0", `18` = "Neu",
                               `19` = "B", `20` = "B", `21` = "Mast2", `22` = "Mac0", `23` = "Neu",
                               `24` = "B", `25` = "B", `26` = "Mast2", `27` = "Mac0", `28` = "Neu",
                               `29` = "B", `30` = "B", `31` = "Mast2", `32` = "Mac0", `33` = "Neu",
                               `34` = "B", `35` = "B")

  # Cell_Type_Order.set <- c("T", "CD4+T", "CD8+T", "B" , "Mac0", "Mac1", "Mac2", "Mac3", "Mast1", "Mast2", "NK", "Neu", "Ery")




##### Save RData #####
  save.image("SeuratObject_CDS_PRJCA001063_ReAnno.RData")




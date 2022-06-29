## SingleRBook Ref: http://bioconductor.org/books/release/SingleRBook/
## Example Ref: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(300000)


##### Load Packages #####
  if(!require("Seurat")) install.packages("Seurat")
  if(!require("tidyverse")) install.packages("tidyverse")
  if(!require("ggpubr")) install.packages("ggpubr")

  library(ggpubr)
  library(tidyverse)
  library(Seurat)

  #### BiocManager installation ####
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("SingleR","scRNAseq","celldex","scran","scater","scuttle")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

##### Function setting #####
  ## Call function
  source("FUN_Anno_SingleR.R")

##### Current path and new folder setting* #####
  ProjectName = paste0("CTAnno_singleR_PRJCA001063S")
  Sampletype = "PDAC"

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

##### Parameter setting* #####
  Remark1 <- "PredbyCTDB" # c("PredbyCTDB","PredbyscRNA")
  RefType <- "BuiltIn_celldex" # c("BuiltIn_celldex","BuiltIn_scRNA")
  celldexDatabase <- "HumanPrimaryCellAtlasData"
  # c("BlueprintEncodeData","DatabaseImmuneCellExpressionData","HumanPrimaryCellAtlasData","ImmGenData",
  #   "MonacoImmuneData","MouseRNAseqData","NovershternHematopoieticData")
  de.method <- "classic"

  ## Parameter of classifySingleR
  quantile = 0.8
  tune.thresh = 0.05
  sd.thresh = 1

  Remark <- paste0(Remark1,"_",de.method,"_",
                   "qua",quantile,"_tun",tune.thresh,"_sd",sd.thresh)

  SmallTest = TRUE

#### Load data #####
  load("SeuratObject_CDS_PRJCA001063.RData")

  ## SeuObj_Ref
  scRNA.SeuObj_Ref <- scRNA.SeuObj

  if(SmallTest == TRUE){
    ## SeuObj_Ref for small test
    # CTFeatures.SeuObj <- scRNA.SeuObj_Ref[,scRNA.SeuObj_Ref$CELL %in% sample(scRNA.SeuObj_Ref$CELL,1000)] ## For small test
    CTFeatures.SeuObj <- scRNA.SeuObj_Ref[,scRNA.SeuObj_Ref@meta.data[[1]] %in% sample(scRNA.SeuObj_Ref@meta.data[[1]],1000)] ## For small test
    ## SeuObj_Tar for small test
    # scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj$CELL %in% sample(scRNA.SeuObj$CELL,1000)] ## For small test
    scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[[1]] %in% sample(scRNA.SeuObj@meta.data[[1]],1000)] ## For small test
  }else{
    ## SeuObj_Ref for full data
    CTFeatures.SeuObj <- scRNA.SeuObj_Ref
  }


##### Run singleR #####
  #### Presetting ####
  SingleRResult.lt <- Anno_SingleR(scRNA.SeuObj, RefType = RefType, celldexDatabase = celldexDatabase,
                                   CTFeatures.SeuObj = CTFeatures.SeuObj,
                                   quantile = quantile, tune.thresh = tune.thresh, sd.thresh = sd.thresh,
                                   de.method = de.method,
                                   Remark = Remark, Save.Path = paste0(Save.Path,"/",Remark), ProjectName = "CT")

  #### Try Parameter ####
  CC_Anno.df <- as.data.frame(matrix(nrow=0, ncol=7))
  colnames(CC_Anno.df) <- c("TestID", "Tool", "Type","Set", "quantile", "tune_Thr","SD_Thr")

  #### PredbyCTDB ####
  for (i in seq(0.6,1,0.2)) {
    for (j in seq(0.03,0.07,0.01)) {
      for (k in c(1,2)) {
        Remark1 <- "PredbyCTDB"
        de.method <- "classic"
        RefType <- "BuiltIn_celldex"
        Remark <- paste0(Remark1,"_",de.method,"_",
                         "qua",i,"_tun",j,"_sd",k)

        SingleRResult.lt <- Anno_SingleR(scRNA.SeuObj, RefType = RefType, celldexDatabase = "HumanPrimaryCellAtlasData",
                                         quantile = i, tune.thresh = j, sd.thresh = k,
                                         Remark = Remark,Save.Path = paste0(Save.Path,"/",Remark), ProjectName = "CT")
        scRNA.SeuObj <- SingleRResult.lt[["scRNA.SeuObj"]]

        CC_Anno_Temp.df <- data.frame(TestID = "Predict", Tool = "singleR", Type = "PDAC",
                                      Set = Remark1, quantile = i, tune_Thr = j, SD_Thr = k)
        CC_Anno.df <- rbind(CC_Anno.df, CC_Anno_Temp.df)
      }
    }
  }
  rm(i,j,k,CC_Anno_Temp.df, Remark, Remark1, de.method, RefType)

    #### Create check dataframe ####
    CC.df <- scRNA.SeuObj@meta.data[,(ncol(scRNA.SeuObj@meta.data)-nrow(CC_Anno.df)+1):ncol(scRNA.SeuObj@meta.data)]
    CC_Anno.df$TestID <- colnames(CC.df)

    # TTT <- gsub(CC.df, pattern = " ", replacement = "_")
    # TTT <- sub(" ", "_", CC.df)
    # TTT <- CC.df
    CC.df <- lapply(CC.df, gsub, pattern = "_", replacement = " ", fixed = TRUE) %>%
             lapply(., gsub, pattern = "cells", replacement = "cell", fixed = TRUE) %>%
             lapply(., gsub, pattern = "Macrophage", replacement = "Macrophage cell", fixed = TRUE) %>%
             lapply(., gsub, pattern = "Fibroblasts", replacement = "Fibroblast cell", fixed = TRUE) %>%
             lapply(., gsub, pattern = "Epithelial cell", replacement = "Ductal cell type 1", fixed = TRUE) %>%
             as.data.frame()

    CC_CT.df <- data.frame(Cell_type = scRNA.SeuObj@meta.data[,"Cell_type"])
    CC.df <- cbind(CC_CT.df, CC.df)
    rm(CC_CT.df)

    CTReplace <- function(CC.df,ColN=1, ReferCT = "Actual",Replacement="Other") {
      CC.df[!CC.df[,ColN] %in% c(CC.df[,ReferCT] %>% unique()),][,ColN] <- Replacement
      return(CC.df)
    }

    for (i in 2:ncol(CC.df)) {
      CC.df <- CTReplace(CC.df,ColN=i, ReferCT = "Cell_type",Replacement="Other")
    }


  #### PredbyscRNA ####
  CC_Anno2.df <- as.data.frame(matrix(nrow=0, ncol=7))
  colnames(CC_Anno2.df) <- c("TestID", "Tool", "Type","Set", "quantile", "tune_Thr","SD_Thr")

  for (i in seq(0.6,1,0.2)) {
    for (j in seq(0.03,0.07,0.01)) {
      for (k in c(1,2)) {
        Remark1 <- "PredbyscRNA"
        de.method <- "classic"
        RefType <- "BuiltIn_scRNA"
        Remark <- paste0(Remark1,"_",de.method,"_",
                         "qua",i,"_tun",j,"_sd",k)

        SingleRResult.lt <- Anno_SingleR(scRNA.SeuObj, RefType = RefType, celldexDatabase = "HumanPrimaryCellAtlasData",
                                         quantile = i, tune.thresh = j, sd.thresh = k,CTFeatures.SeuObj = CTFeatures.SeuObj,
                                         Remark = Remark, Save.Path = paste0(Save.Path,"/",Remark), ProjectName = "CT")
        scRNA.SeuObj <- SingleRResult.lt[["scRNA.SeuObj"]]

        CC_Anno_Temp.df <- data.frame(TestID = "Predict", Tool = "singleR", Type = "PDAC",
                                      Set = Remark1, quantile = i, tune_Thr = j, SD_Thr = k)
        CC_Anno2.df <- rbind(CC_Anno2.df, CC_Anno_Temp.df)
      }
    }
  }
  rm(i,j,k,CC_Anno_Temp.df, Remark, Remark1, de.method, RefType)

    #### Create check dataframe ####
    CC2.df <- scRNA.SeuObj@meta.data[,(ncol(scRNA.SeuObj@meta.data)-nrow(CC_Anno2.df)+1):ncol(scRNA.SeuObj@meta.data)]
    CC_Anno2.df$TestID <- colnames(CC2.df)
    CC_Anno.df <- rbind(CC_Anno.df,CC_Anno2.df)

    CC_CT.df <- data.frame(Cell_type = scRNA.SeuObj@meta.data[,"Cell_type"])
    CC.df <- cbind(CC.df, CC2.df)
    rm(CC_CT.df,CC2.df,CC_Anno2.df)

    CC.df$Cell_type <- as.character(CC.df$Cell_type)
    CC.df <- rbind(CC.df,"Other")

##### Verification (CellCheck) #####
  #### Install ####
  ## Check whether the installation of those packages is required
  Package.set <- c("tidyverse","caret","cvms","DescTools","devtools","ggthemes")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  # library(Seurat)
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  ## install CellCheck
  # Install the CellCheck package
  detach("package:CellCheck", unload = TRUE)
  devtools::install_github("Charlene717/CellCheck")
  # Load CellCheck
  library(CellCheck)

  #### Run CellCheck ####

  ## For one prediction
  DisCMSet.lt = list(Mode = "One", Actual = "Cell_type", Predict = "singleR_PredbyCTDB_classic_qua0.6_tun0.03_sd1" , FilterSet1 = "Tool", FilterSet2 = "singleR" , Remark = "") # Mode = c("One","Multiple")
  BarChartSet.lt <- list(Mode = "One", Metrics = "Balanced.Accuracy", XValue = "Set", Group = "Tool", Remark = "")
  LinePlotSet.lt <- list(Mode = "One", Metrics = "Balanced.Accuracy", XValue = "tune_Thr", Group = "Set", Remark = "")
  CCR_cm_DisMult.lt <- CellCheck_DisMult(CC.df, CC_Anno.df,
                                         DisCMSet.lt = DisCMSet.lt,
                                         BarChartSet.lt = BarChartSet.lt,
                                         LinePlotSet.lt = LinePlotSet.lt,
                                         Save.Path = Save.Path, ProjectName = paste0("CellCheck_",ProjectName))

  ## For multiple prediction
  DisCMSet.lt = list(Mode = "Multiple", Actual = "Cell_type", FilterSet1 = "Tool", FilterSet2 = "singleR" , Remark = "_All") # Mode = c("One","Multiple")
  BarChartSet.lt <- list(Mode = "Multiple", XValue = "Set", Group = "Tool", Remark = "_All")
  LinePlotSet.lt <- list(Mode = "Multiple", XValue = "tune_Thr", Group = "Set", Remark = "_All")
  Sum_DisMult.df <- CellCheck_DisMult(CC.df, CC_Anno.df,
                                      DisCMSet.lt = DisCMSet.lt,
                                      BarChartSet.lt = BarChartSet.lt,
                                      LinePlotSet.lt=LinePlotSet.lt,
                                      Save.Path = Save.Path, ProjectName = paste0("CellCheck_",ProjectName))



##### Session information #####
  sessionInfo()
  ## Ref: https://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
  writeLines(capture.output(sessionInfo()), paste0(Save.Path,"/sessionInfo.txt"))

##### Save RData #####
  save.image(paste0(Save.Path,"/SeuratObject_",ProjectName,".RData"))





# #### Old Version ####
#   #### Run SingleR ####
#   SingleRResult.lt <- Anno_SingleR(scRNA.SeuObj, RefType = "BuiltIn_celldex", celldexDatabase = "HumanPrimaryCellAtlasData",
#                                    quantile = quantile, tune.thresh = tune.thresh, sd.thresh = sd.thresh,
#                                    Remark = "PredbyCTDB",Save.Path = Save.Path, ProjectName = ProjectName)
#
#   scRNA.SeuObj <- SingleRResult.lt[["scRNA.SeuObj"]]
#   SingleRResult2.lt <- Anno_SingleR(scRNA.SeuObj, RefType = "BuiltIn_scRNA", celldexDatabase = "HumanPrimaryCellAtlasData",
#                                    quantile = quantile, tune.thresh = tune.thresh, sd.thresh = sd.thresh,
#                                    Remark = "PredbyscRNA",CTFeatures.SeuObj = CTFeatures.SeuObj, de.method = "classic",
#                                    Save.Path = Save.Path, ProjectName = ProjectName)
#   scRNA.SeuObj <- SingleRResult2.lt[["scRNA.SeuObj"]]
#
#
#
# ##### Verification (CellCheck) #####
#   #### Install ####
#   ## Check whether the installation of those packages is required
#   Package.set <- c("tidyverse","caret","cvms","DescTools","devtools","ggthemes")
#   for (i in 1:length(Package.set)) {
#     if (!requireNamespace(Package.set[i], quietly = TRUE)){
#       install.packages(Package.set[i])
#     }
#   }
#   ## Load Packages
#   # library(Seurat)
#   lapply(Package.set, library, character.only = TRUE)
#   rm(Package.set,i)
#
#   ## install CellCheck
#   # Install the CellCheck package
#   detach("package:CellCheck", unload = TRUE)
#   devtools::install_github("Charlene717/CellCheck")
#   # Load CellCheck
#   library(CellCheck)
#
#   #### Run CellCheck ####
#   ## Create check dataframe
#   CC.df <- scRNA.SeuObj@meta.data[,c("Cell_type","singleR_classic_PredbyscRNA", "singleR_classic_PredbyCTDB")]
#
#   CC.df <- data.frame(lapply(CC.df, as.character), stringsAsFactors=FALSE)
#
#   colnames(CC.df) <- c("Actual","Predict1","Predict2")
#   #CC.df$Actual <- as.character(CC.df$Actual)
#
#   CC.df$Predict2 <- gsub("_", " ", CC.df$Predict2)
#   CC.df$Predict2 <- gsub("cells", "cell", CC.df$Predict2)
#   CC.df$Predict2 <- gsub("Macrophage", "Macrophage cell", CC.df$Predict2)
#   CC.df$Predict2 <- gsub("Fibroblasts", "Fibroblast cell", CC.df$Predict2)
#   CC.df$Predict2 <- gsub("Epithelial cell", "Ductal cell type 1", CC.df$Predict2)
#
#   CC.df[!CC.df$Predict2 %in% c(CC.df$Actual %>% unique()),]$Predict2 <- "Other"
#   # CC.df <- rbind(CC.df,"NotMatch")  #CC.df[nrow(CC.df)+1,1:ncol(CC.df)] <- "Other"
#
#   CC_Anno.df <- data.frame(TestID = c("Predict1","Predict2"),
#                            Tool = "singleR",
#                            Type = "PDAC",
#                            Set = c("singleRPredbyscRNA", "singleRPredbyCTDB"))
#
#   ## For one prediction
#   ## For one prediction
#   DisCMSet.lt = list(Mode = "One", Actual = "Actual", Predict = "Predict1" , FilterSet1 = "Tool", FilterSet2 = "singleR" , Remark = "") # Mode = c("One","Multiple")
#   BarChartSet.lt <- list(Mode = "One", Metrics = "Balanced.Accuracy", XValue = "Set", Group = "Tool", Remark = "")
#   LinePlotSet.lt <- list(Mode = "One", Metrics = "Balanced.Accuracy", XValue = "Set", Group = "Tool", Remark = "")
#   CCR_cm_DisMult.lt <- CellCheck_DisMult(CC.df, CC_Anno.df,
#                                          DisCMSet.lt = DisCMSet.lt,
#                                          BarChartSet.lt = BarChartSet.lt,
#                                          LinePlotSet.lt = LinePlotSet.lt,
#                                          Save.Path = Save.Path, ProjectName = paste0("CellCheck_",ProjectName))
#
#
#
# ##### Save RData #####
#   save.image(paste0(Save.Path,"/SeuratObject_",ProjectName,".RData"))
#
#
#
#

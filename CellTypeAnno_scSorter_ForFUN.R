## Ref: https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages #####
  #### Basic installation ####
  Package.set <- c("tidyverse","scSorter","Seurat","stringr","ggpubr","magrittr","dplyr")
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
  source("FUN_Anno_scSorter.R")

##### Current path and new folder setting* #####
  ProjectName = paste0("CTAnno_scSorter_PRJCA001063S")
  Sampletype = "PDAC"

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

##### Parameter setting* #####
  Remark1 <- "PredbyscRNA" # c("PredbyCTDB","PredbyscRNA")
  RefType <- "BuiltIn_scRNA" # c("BuiltIn_celldex","BuiltIn_scRNA")
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


##### Set References #####
  ## Create cell type markers dataframe
  Idents(scRNA.SeuObj) <- scRNA.SeuObj$Cell_type
  CellType.markers.df <- FindAllMarkers(scRNA.SeuObj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)



  scRNA.SeuObj <- Anno_scSorter(scRNA.SeuObj, CellType.markers.df,
                                Path = Save.Path, projectName = ProjectName,
                                log2FC_CN = "avg_log2FC",log2FC_Thr = 1 ,
                                RefName = "Cell_type",
                                pValue_CN = "p_val", pVal_Thr = 0.05,
                                Gene_CN = "gene",
                                Cluster_CN = "cluster")

  ## Create anno.df
  CTTop.markers <-  CellType.markers.df %>%
                    group_by(cluster) %>%
                    top_n(n = 10, wt = avg_log2FC)

  #DoHeatmap(scRNA.SeuObj, features = CTTop.markers$gene) + NoLegend()
  log2FC_CN <- "avg_log2FC"
  Gene_CN <- "gene"
  Cluster_CN <- "cluster"

  anno.df <- data.frame(Type = CellType.markers.df[,Cluster_CN],
                        Marker = CellType.markers.df[,Gene_CN],
                        Weight = CellType.markers.df[,log2FC_CN])

##### Set Target Gene expression matrix #####
  GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") # normalized data matrix
  # GeneExp.df <- scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame()

#### Run scSorter ####
  library(scSorter)

  scSorter.obj <- scSorter(GeneExp.df,anno.df)

  scRNA.SeuObj$scSorterPred <- scSorter.obj[["Pred_Type"]]
  p.CTPred1 <- DimPlot(scRNA.SeuObj, reduction = "umap", group.by ="scSorterPred" ,label = TRUE, pt.size = 0.5) + NoLegend()
  p.CTPred1
  p.CT1 <- DimPlot(scRNA.SeuObj, reduction = "umap", group.by ="Cell_type" ,label = TRUE, pt.size = 0.5) + NoLegend()
  p.CT1

  library(ggpubr)
  p.CTComp1 <- ggarrange(p.CT1, p.CTPred1, common.legend = TRUE, legend = "top")
  p.CTComp1

  scRNA.SeuObj@meta.data[[paste0("scSorterPred","_LogFC",1,"_pV",0.05)]] <- scSorter.obj[["Pred_Type"]]
  DimPlot(scRNA.SeuObj, reduction = "umap",
          group.by = paste0("scSorterPred","_LogFC",1,"_pV",0.05),label = TRUE, pt.size = 0.5) + NoLegend()
  DimPlot(scRNA.SeuObj, reduction = "umap", group.by ="celltype" ,label = TRUE, pt.size = 0.5) + NoLegend()


##### Save RData #####
  save.image(paste0(Save.Path,"/SeuratObject_",ProjectName,".RData"))





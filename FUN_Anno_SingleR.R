## SingleRBook Ref: http://bioconductor.org/books/release/SingleRBook/
## Example Ref: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html


Anno_SingleR <- function(scRNA.SeuObj, RefType = "BuiltIn_celldex", celldexDatabase = "HumanPrimaryCellAtlasData",
                         Remark = "PredbyCTDB",Save.Path = "", ProjectName = "",
                         RefName = "Cell_type",
                         quantile = 0.8, tune.thresh = 0.05, sd.thresh = 1,
                         CTFeatures.SeuObj = CTFeatures.SeuObj, de.method = "classic"
                         ) {
## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

##### Set References #####
  if(RefType == "BuiltIn_celldex"){
    #### Database: Bulk reference setting for Cell type features ####
    library(celldex)

    if(celldexDatabase == "BlueprintEncodeData"){
      CTFeatures <- BlueprintEncodeData()
    }else if(celldexDatabase == "DatabaseImmuneCellExpressionData"){
      CTFeatures <- DatabaseImmuneCellExpressionData()
    }else if(celldexDatabase == "HumanPrimaryCellAtlasData"){
      CTFeatures <- HumanPrimaryCellAtlasData()
    }else if(celldexDatabase == "ImmGenData"){
      CTFeatures <- ImmGenData()
    }else if(celldexDatabase == "MonacoImmuneData"){
      CTFeatures <- MonacoImmuneData()
    }else if(celldexDatabase == "MouseRNAseqData"){
      CTFeatures <- MouseRNAseqData()
    }else if(celldexDatabase == "NovershternHematopoieticData"){
      CTFeatures <- NovershternHematopoieticData()
    }else{
      print("Error in database setting!")
    }

    #### Demo dataset ####
    # library(celldex)
    # hpca.se <- HumanPrimaryCellAtlasData()
    # hpca.se
  }else if(RefType =="BuiltIn_scRNA"){
    #### single-cell reference setting for Cell type features ####
    ## Prepossessing
    CTFeatures <- as.SingleCellExperiment(CTFeatures.SeuObj)
    CTFeatures$label <- CTFeatures@colData@listData[[RefName]]
    CTFeatures <- CTFeatures[,!is.na(CTFeatures$label)]
    # CTFeatures <- logNormCounts(CTFeatures)
    #rm(CTFeatures.SeuObj)

    #### Demo dataset ####
    # library(scRNAseq)
    # sceM <- MuraroPancreasData()
    #
    # # One should normally do cell-based quality control at this point, but for
    # # brevity's sake, we will just remove the unlabelled libraries here.
    # sceM <- sceM[,!is.na(sceM$label)]
    #
    # # SingleR() expects reference datasets to be normalized and log-transformed.
    # library(scuttle)
    # sceM <- logNormCounts(sceM)

  }

##### Set Target SeuObj #####
  ## Prepossessing
  scRNA <- as.SingleCellExperiment(scRNA.SeuObj)
  #### Demo dataset ####
  # library(scRNAseq)
  # hESCs <- LaMannoBrainData('human-es')
  # hESCs <- hESCs[,colSums(counts(hESCs)) > 0] # Remove libraries with no counts.
  # hESCs <- logNormCounts(hESCs)
  # hESCs <- hESCs[,1:100]

#### Run SingleR ####
  library(SingleR)
  if(RefType == "BuiltIn_celldex"){
    SingleR.lt <- SingleR(test = scRNA, ref = CTFeatures, assay.type.test=1,
                          quantile = quantile, tune.thresh = tune.thresh, sd.thresh = sd.thresh,
                          labels = CTFeatures$label.main , de.method= de.method)#, de.method="wilcox") #  de.method = c("classic", "wilcox", "t")

  }else if(RefType == "BuiltIn_scRNA"){
    SingleR.lt <- SingleR(test = scRNA, ref = CTFeatures, assay.type.test=1,
                          quantile = quantile, tune.thresh = tune.thresh, sd.thresh = sd.thresh,
                          labels = CTFeatures$label , de.method= de.method)#, de.method="wilcox") #  de.method = c("classic", "wilcox", "t")
  }

  SingleR.lt

  # Summarizing the distribution:
  CTCount_Sum.df <- table(SingleR.lt$labels) %>%
                       as.data.frame() %>%
                       dplyr::rename(Cell_Type = Var1, Count = Freq)

##### Annotation diagnostics #####
  ## Overall
  p.ScoreHeatmap1 <- plotScoreHeatmap(SingleR.lt)
  p.ScoreHeatmap1
  p.DeltaDist1 <- plotDeltaDistribution(SingleR.lt, ncol = 5)
  p.DeltaDist1
  summary(is.na(SingleR.lt$pruned.labels))

  pdf(file = paste0(Save.Path,"/",ProjectName,"_",Remark,"_AnnoDiag.pdf"),
      width = 10,  height = 7
  )
  p.ScoreHeatmap1 %>% print()
  p.DeltaDist1 %>% print()
  dev.off()



  all.markers <- metadata(SingleR.lt)$de.genes
  scRNA@colData@listData[[paste0("labels_",Remark)]] <- SingleR.lt$labels ## scRNA$labels <- SingleR.lt$labels


  ## HeatmapCTmarkers
  # # Endothelial cell-related markers
  # library(scater)
  # plotHeatmap(scRNA, order_columns_by="labels",
  #             features = unique(unlist(all.markers[["Endothelial_cells"]])))

  # ## All in one PDF
  # library(scater)
  # pdf(file = paste0(Save.Path,"/",ProjectName,"_",Remark,"_HeatmapCTmarkers.pdf"),
  #     width = 12,  height = 7
  # )
  #
  #   for (i in 1:length(all.markers)) {
  #     p <-  plotHeatmap(scRNA, order_columns_by = paste0("labels_",Remark),
  #                       features=unique(unlist(all.markers[[i]])))
  #     print(p)
  #   }
  #
  # dev.off() # graphics.off()

  ## Split PDF file
  library(scater)
  for (i in 1:length(all.markers)) {
  pdf(file = paste0(Save.Path,"/",ProjectName,"_",Remark,"_HeatmapCTmarkers_",names(all.markers[i]),".pdf"),
      width = 12,  height = 7
  )


    p <-  plotHeatmap(scRNA, order_columns_by = paste0("labels_",Remark),
                      features=unique(unlist(all.markers[[i]])))
    print(p)

  dev.off() # graphics.off()
  }

  ## Plot UMAP
  scRNA.SeuObj@meta.data[[paste0("singleR_",Remark)]]<- SingleR.lt$labels # scRNA.SeuObj$singleRPredbyCTDB <- SingleR.lt$labels
  p.CTPred1 <- DimPlot(scRNA.SeuObj, reduction = "umap", group.by = paste0("singleR_",Remark) ,label = TRUE, pt.size = 0.5) # + NoLegend()
  p.CTPred1
  p.CT1 <- DimPlot(scRNA.SeuObj, reduction = "umap", group.by = RefName ,label = TRUE, pt.size = 0.5) # + NoLegend()
  p.CT1

  library(ggpubr)
  p.CTComp1 <- ggarrange(p.CT1, p.CTPred1, common.legend = TRUE, legend = "top")
  p.CTComp1


  pdf(file = paste0(Save.Path,"/",ProjectName,"_",Remark,"_CompareCTUMAP.pdf"),
      width = 12,  height = 7
  )
    p.CTComp1 %>% print()
    p.CTPred1 %>% print()
    p.CT1 %>% print()
  dev.off()


  Output.lt <- list(SingleR.lt = SingleR.lt,
                    scRNA.SeuObj = scRNA.SeuObj,
                    SingleCellExperiment = scRNA,
                    Summarize.df = CTCount_Sum.df
                    )

  return(Output.lt)


}


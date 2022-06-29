## Ref: https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html
Anno_scSorter <- function(scRNA.SeuObj, CTFilter.Markers.df,
                          Path = "", projectName = "",
                          RefName = "Cell_type",
                          log2FC_CN = "avg_log2FC",log2FC_Thr = 1 ,
                          pValue_CN = "p_val", pVal_Thr = 0.05,
                          Gene_CN = "gene",
                          Cluster_CN = "cluster") {

  ##### Section 1 - Preliminaries #####
  ## Presetting
  memory.limit(300000)

  ## Instal and Load Packages
    ## Basic installation
    Package.set <- c("tidyverse","scSorter","Seurat","stringr","magrittr","dplyr")
    ## Check whether the installation of those packages is required from basic
    for (i in 1:length(Package.set)) {
      if (!requireNamespace(Package.set[i], quietly = TRUE)){
        install.packages(Package.set[i])
      }
    }
    ## Load Packages
    lapply(Package.set, library, character.only = TRUE)
    rm(Package.set,i)



  ##### Section 2 - Preprocessing the data #####

    ## Filter cell marker
    CTFilter.Markers.df <-  CTFilter.Markers.df %>%
                            dplyr::filter(.,.[,log2FC_CN] >= log2FC_Thr) %>%
                            dplyr::filter(.,.[pValue_CN] <= pVal_Thr)
    ## print Heatmap
    p.Heatmap <- DoHeatmap(scRNA.SeuObj, features = CTFilter.Markers.df$gene,
                           angle = 90, size = 3) +
                  scale_fill_gradient2(low = "#5283ff",
                                       mid = "white",
                                       high = "#ff5c5c") +
                  theme(axis.text.y = element_text(size  = 5)) +
                  theme(legend.position = "bottom" )
    print(p.Heatmap)

    ## Export pdf
    pdf(
      file = paste0(Path, "/",projectName,"_Heatmap_","_LogFC",log2FC_Thr,"_pV",pVal_Thr,".pdf"),
      width = 10,  height = 8
    )

     print(p.Heatmap)

    dev.off()
    ## Create anno.df
    anno.df <- data.frame(Type = CTFilter.Markers.df[,Cluster_CN],
                          Marker = CTFilter.Markers.df[,Gene_CN],
                          Weight = CTFilter.Markers.df[,log2FC_CN])

    ## Export Gene expression matrix
    GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") # normalized data matrix

  ##### Section 3 - Running scSorter #####
    scSorter.obj <- scSorter(GeneExp.df, anno.df)
    scRNA.SeuObj@meta.data[[paste0("scSorterPred","_LogFC",log2FC_Thr,"_pV",pVal_Thr)]] <- scSorter.obj[["Pred_Type"]]

    ## print UMAP
    ## Export pdf
    pdf(
      file = paste0(Path, "/",projectName,"_UMAP_","_LogFC",log2FC_Thr,"_pV",pVal_Thr,".pdf"),
      width = 7,  height = 7
    )
      p.UMAP1 <- DimPlot(scRNA.SeuObj, reduction = "umap",
                         group.by = paste0("scSorterPred","_LogFC",log2FC_Thr,"_pV",pVal_Thr),label = TRUE, pt.size = 0.5) + NoLegend()
      print(p.UMAP1)

      p.UMAP2 <- DimPlot(scRNA.SeuObj, reduction = "umap",
                         group.by = RefName ,label = TRUE, pt.size = 0.5) + NoLegend()
      print(p.UMAP2)
    dev.off()

    return(scRNA.SeuObj)

}

Extract_Feature <- function(CellType.markers.df,
                            TOPN = "",
                            log2FC_CN = "avg_log2FC",
                            log2FC_Thr = 1,
                            pValue_CN = "p_val",
                            pValue_Thr = 0.05,
                            pValueAdj_CN = "p_val_adj",
                            pValueAdj_Thr = 0.05

                            ) {
  ##### Presetting ######
    memory.limit(300000)

  ##### Load Packages #####
    #### Basic installation ####
    Package.set <- c("tidyverse", "stringr", "dplyr")
    ## Check whether the installation of those packages is required from basic
    for (i in 1:length(Package.set)) {
      if (!requireNamespace(Package.set[i], quietly = TRUE)){
        install.packages(Package.set[i])
      }
    }
    ## Load Packages
    lapply(Package.set, library, character.only = TRUE)
    rm(Package.set,i)

  #### Filter marker ####
    if(TOPN != ""){
      ## By Rank
      CTFilter.Markers.df <-  CellType.markers.df %>%
                              # group_by(.,.[,Cluster]) %>% ## Bug
                              # top_n(.,n = TOPN, wt = .[,log2FC_CN]) ## Error
                              group_by(cluster) %>%
                              top_n(.,n = TOPN, wt = avg_log2FC) ## Bug

    }else{
      ## By filters
      CTFilter.Markers.df <-  CellType.markers.df %>%
                              dplyr::filter(.,.[,log2FC_CN] >= log2FC_Thr) %>%
                              dplyr::filter(.,.[pValue_CN] <= pValue_Thr) %>%
                              dplyr::filter(.,.[,pValueAdj_CN] <= pValueAdj_Thr)

    }
    return(CTFilter.Markers.df)

}



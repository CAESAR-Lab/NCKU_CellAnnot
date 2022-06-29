## SingleRBook Ref: http://bioconductor.org/books/release/SingleRBook/
## Example Ref: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)


##### Load Packages #####
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


##### Using built-in references #####
  library(celldex)
  hpca.se <- HumanPrimaryCellAtlasData()
  hpca.se

  library(scRNAseq)
  hESCs <- LaMannoBrainData('human-es')
  hESCs <- hESCs[,1:100]

  library(SingleR)
  pred.hesc <- SingleR(test = hESCs, ref = hpca.se, assay.type.test=1,
                       labels = hpca.se$label.main)

  pred.hesc

  # Summarizing the distribution:
  table(pred.hesc$labels)

##### Using single-cell references   #####
  library(scRNAseq)
  sceM <- MuraroPancreasData()

  # One should normally do cell-based quality control at this point, but for
  # brevity's sake, we will just remove the unlabelled libraries here.
  sceM <- sceM[,!is.na(sceM$label)]

  # SingleR() expects reference datasets to be normalized and log-transformed.
  library(scuttle)
  sceM <- logNormCounts(sceM)

  sceG <- GrunPancreasData()
  sceG <- sceG[,colSums(counts(sceG)) > 0] # Remove libraries with no counts.
  sceG <- logNormCounts(sceG)
  sceG <- sceG[,1:100]

  library(SingleR)
  pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, de.method="wilcox")
  table(pred.grun$labels)

##### Annotation diagnostics #####
  plotScoreHeatmap(pred.grun)
  plotDeltaDistribution(pred.grun, ncol = 3)
  summary(is.na(pred.grun$pruned.labels))

  all.markers <- metadata(pred.grun)$de.genes
  sceG$labels <- pred.grun$labels

  # Beta cell-related markers
  library(scater)
  plotHeatmap(sceG, order_columns_by="labels",
              features=unique(unlist(all.markers$beta)))


##### Session information #####
  sessionInfo()



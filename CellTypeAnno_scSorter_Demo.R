## Ref: https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html

# ##### Presetting ######
#   rm(list = ls()) # Clean variable
#   memory.limit(150000)


##### Section 1 - Preliminaries #####
  library(scSorter)

  ## Examing the data
  load(url('https://github.com/hyguo2/scSorter/blob/master/inst/extdata/TMpancreas.RData?raw=true'))
  expr[1:5, 1:5]
  dim(expr)

##### Section 2 - Preprocessing the data #####

  topgenes = xfindvariable_genes(expr, ngenes = 2000)

  expr = xnormalize_scData(expr)

  topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
  topgenes = topgenes[topgene_filter]


  picked_genes = unique(c(anno$Marker, topgenes))
  expr = expr[rownames(expr) %in% picked_genes, ]

  ## Create anno.df
  CellType.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> CTTop.markers
  #DoHeatmap(scRNA.SeuObj, features = CTTop.markers$gene) + NoLegend()
  anno.df <- data.frame(Type = CTTop.markers$cluster,
                        Marker = CTTop.markers$gene,
                        Weight = CTTop.markers$avg_log2FC)



##### Section 3 - Running scSorter #####
  rts <- scSorter(expr, anno)

  # Viewing Results
  print(table(rts$Pred_Type))

  mis_rate = 1 - mean(rts$Pred_Type == true_type)
  round(mis_rate, 4)

  table(true_type, rts$Pred_Type)


# References
# Tabula Muris Consortium and others. 2018. “Single-Cell Transcriptomics of 20 Mouse Organs Creates a Tabula Muris.” Nature 562 (October). Genome Research: 367–72. doi:10.1038/s41586-018-0590-4.

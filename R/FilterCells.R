#' Filter Seurat Object automatically
#'
#' The number of features and UMIs (nFeature_RNA and nCount_RNA) are automatically calculated for every object by Seurat.
#' For non-UMI data, nCount_RNA represents the sum of the non-normalized values within a cell
#' We calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`.
#' We use raw count data since this represents non-transformed and non-log-normalized counts
#' The % of UMI mapping to MT-features is a common scRNA-seq QC metric.
#'
#' @author David John & Mariano Ruz Jurado
#' @param seuratObject a Seurat Object
#' @return filtered seurat object
#' @export
FilterDeadCellsByQuantile <- function(seuratObject, lowQuantile=0.1 , highQuantile=0.95, maxMito=0.1){
  sizeBefore<-length(seuratObject@meta.data$orig.ident)
  cat("FilterByQuantile\n")
  #For some unknown reasons these variables need to be global for the subset function, otherwise there is an eval() unknown variable error
  lowQuantile<<-lowQuantile
  highQuantile<<-highQuantile
  maxMito<<-maxMito
  sample<-unique(seuratObject$sample)
  Quality <- data.frame(UMI=seuratObject$nCount_RNA, nGene=seuratObject$nFeature_RNA, label = factor(seuratObject$sample), percent.mito=seuratObject$percent.mito)

  Quantile.low.UMI <- Quality |> dplyr::group_by(label) |>
    dplyr::summarise(UMI = list(enframe(quantile(UMI,probs = lowQuantile)))) |>
    tidyr::unnest(cols = c(UMI))

  Quantile.high.UMI <- Quality |> dplyr::group_by(label) |>
    dplyr::summarise(UMI = list(enframe(quantile(UMI,probs = highQuantile)))) |>
    tidyr::unnest(cols = c(UMI))

  Quantile.low.Gene <- Quality |> dplyr::group_by(label) |>
    dplyr::summarise(nGene = list(enframe(quantile(nGene,probs = lowQuantile)))) |>
    tidyr::unnest(cols = c(nGene))

  Quantile.high.Gene <- Quality |> dplyr::group_by(label) |>
    dplyr::summarise(nGene = list(enframe(quantile(nGene,probs = highQuantile)))) |>
    tidyr::unnest(cols = c(nGene))


  df<-seuratObject@meta.data

  gg1<- ggplot2::ggplot(Quality, aes(x="nUMI", y=UMI)) + geom_violin(scale = "width") +
    theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none", axis.text.x = element_text(size=12, face = "bold"), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) +
    geom_hline(yintercept = Quantile.high.UMI$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.high.UMI$value, label=Quantile.high.UMI$value , vjust = -1)) +
    geom_hline(yintercept = Quantile.low.UMI$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.low.UMI$value, label=Quantile.low.UMI$value , vjust = -1))

  gg2<- ggplot2::ggplot(Quality, aes(x="nFeature_RNA", y=nGene)) + geom_violin(scale = "width") +
    theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none", axis.text.x = element_text(size=12, face = "bold"), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) +
    geom_hline(yintercept = Quantile.high.Gene$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.high.Gene$value, label=Quantile.high.Gene$value , vjust = -1)) +   geom_hline(yintercept = Quantile.low.Gene$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.low.Gene$value, label=Quantile.low.Gene$value , vjust = -1))


  gg3<- ggplot2::ggplot(Quality, aes(x=" % Mt Content", y=percent.mito)) + geom_violin(scale = "width") +
    theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none", axis.text.x = element_text(size=12, face = "bold"), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) +
    geom_hline(yintercept = maxMito, color="red", linetype="dashed") + geom_text(aes(0.9,maxMito, label=maxMito , vjust = -1))

  gg<-ggpubr::ggarrange(gg1,gg2,gg3, ncol = 3)


  gg<-ggpubr::annotate_figure(gg, fig.lab = sample, fig.lab.pos = "top", fig.lab.size = 15, fig.lab.face = "bold")

  seuratObject<- subset(x= seuratObject, subset = nCount_RNA < Quantile.high.UMI$value & nCount_RNA > Quantile.low.UMI$value &
                          nFeature_RNA < Quantile.high.Gene$value & nFeature_RNA > Quantile.low.Gene$value & percent.mito < maxMito)



  diff<-  sizeBefore -length(seuratObject@meta.data$orig.ident)
  cat("Filtered ",diff, "from" , sizeBefore, " cells\n", "(minFeatures=",Quantile.low.Gene$value, "; maxFeatures=", Quantile.high.Gene$value, "; maxMito=" ,maxMito, ") for ", unique(seuratObject$sample), "\n" )
  rm(maxMito)
  return(list(seuratObject, gg))
}

#' Filter Seurat Object by setting values manually
#'
#' The number of features and UMIs (nFeature_RNA and nCount_RNA) are calculated for every object by Seurat.
#' For non-UMI data, nCount_RNA represents the sum of the non-normalized values within a cell
#' you can set the maximal percentage of mitochondrial features here and store it in object metadata as `percent.mito`.
#' We use raw count data since this represents non-transformed and non-log-normalized counts
#' The % of UMI mapping to MT-features is a common scRNA-seq QC metric.
#'
#' @author David John & Mariano Ruz Jurado
#' @param seuratObject a Seurat object
#' @param minFeatures Number of minimal Features in a cell
#' @param maxFeatures Number of maximal Features in a cell
#' @param minCounts Number of minimal Counts in a cell
#' @param maxCounts Number of maximal Counts in a cell
#' @param maxMito Percentage of maximal mitochondrial content in a cell
#' @return filtered seurat object
#' @export
FilterDeadCells <- function(seuratObject, minFeatures=300, maxFeatures=6000,minCounts=500,maxCounts=15000, maxMito=0.05){
  # The number of features and UMIs (nFeature_RNA and nCount_RNA) are automatically calculated for every object by Seurat.
  # For non-UMI data, nCount_RNA represents the sum of the non-normalized values within a cell
  # We calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`.
  # We use raw count data since this represents non-transformed and non-log-normalized counts
  # The % of UMI mapping to MT-features is a common scRNA-seq QC metric.
  sizeBefore<-length(seuratObject@meta.data$orig.ident)

  #For some unknown reasons these variables need to be global for the subset function, otherwise there is an eval() unknown variable error
  minFeatures<<-minFeatures
  maxFeatures<<-maxFeatures
  maxMito<<-maxMito
  seuratObject <- subset(x = seuratObject, subset = nFeature_RNA > minFeatures & nFeature_RNA < maxFeatures & nCount_RNA > minCounts & nCount_RNA < maxCounts & percent.mito < maxMito)

  diff<-  sizeBefore -length(seuratObject@meta.data$orig.ident)
  percent <- round(diff/sizeBefore*100,digits = 2)
  cat("Filtered ",diff, "from" , sizeBefore, " cells [",percent,"%]\n", "(minFeatures=",minFeatures, "; maxFeatures=", maxFeatures, "; minCounts=" ,minCounts,  "; maxCounts=" ,maxCounts , "; maxMito=" ,maxMito, ") for ", unique(seuratObject$sample), "\n" )
  rm(minFeatures,maxFeatures,maxMito)
  return(seuratObject)
}

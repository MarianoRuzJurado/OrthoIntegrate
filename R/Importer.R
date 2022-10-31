#' Import Single cell sequencing experiments into Seurat and perform normalisation and scale Data and do a summary of mapping stats
#' @author David John & Mariano Ruz Jurado
#' @param pathways A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids Vector of strings that are assigned to the concordant cells
#' @param TenX Logical for TenX result analysis
#' @param performNormalisation Logical if Seurat Normalization should be performed
#' @param performScaling Logical if Seurat ScaleData should be performed
#' @param performVariableGeneDetection Logical if Seurat FindVariableFeatures should be performed
#' @param FilterCells Logical if Cells should be filtered
#' @param FilterByAbsoluteValues Logical if Cells should be filtered by manual values, further arguments must be provided for minFeatures,maxFeatures,minCounts,maxCounts,maxMito
#' @return Merged seurat object
#' @export
Importer <- function(pathway,id, TenX=TRUE, performNormalisation=TRUE, performScaling = FALSE,performVariableGeneDetection=TRUE, FilterCells=TRUE, FilterByAbsoluteValues=FALSE,...) {

  if (TenX) {
    Matrix <- Seurat::Read10X(pathway)
  }  else{
    Matrix <- read.table(pathway,header = TRUE,sep = ",", dec = ".", row.names = 1)
  }
  #catch optional parameters
  optionalParameters <- list(...)

  seuratObject =Seurat::CreateSeuratObject(counts = Matrix, project = id, min.cells = 5)
  seuratObject$sample <- id
  tmp<-unlist(strsplit(id,split = "-"))
  seuratObject$condition <- paste0(tmp[1:length(tmp)-1],collapse = "-")

  mito.features <- grep(pattern = "^MT-", x = rownames(x = seuratObject), value = TRUE)
  if (length(mito.features)<10) {
    mito.features <- grep(pattern = "^mt-", x = rownames(x = seuratObject), value = TRUE)
  }
  if (length(mito.features)<1) {
    warning("Error: Could not find MT genes")

  }

  percent.mito <- Matrix::colSums(x = GetAssayData(object = seuratObject, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObject, slot = 'counts'))
  seuratObject$percent.mito <- percent.mito

  #write QC to file
  p1<-VlnPlot(object = seuratObject, features = c("nFeature_RNA"), ncol = 1, pt.size = 0) + theme(legend.position = "None", axis.title.x = element_blank(), axis.text.x = element_blank())
  p2<-VlnPlot(object = seuratObject, features = c("nCount_RNA"), ncol = 1, pt.size = 0) + theme(legend.position = "None", axis.title.x = element_blank(), axis.text.x = element_blank())
  p3<-VlnPlot(object = seuratObject, features = c("percent.mito"), ncol = 1, pt.size = 0) + theme(legend.position = "None", axis.title.x = element_blank(), axis.text.x = element_blank())
  gg_preFiltering <-  ggpubr::ggarrange(p1,p2,p3, nrow = 1)
  gg_preFiltering <- ggpubr::annotate_figure(gg_preFiltering, top = text_grob(id,face="bold",color = "darkred",size=18,hjust = 0.2))
  ggsave(filename = paste0(pathway,"QC_preFiltered.svg"),device = "svg", width = 10,height = 10)

  if (FilterCells==TRUE) {
    message("start Filtering")
    if (FilterByAbsoluteValues==TRUE) {
      if (is.null(optionalParameters$minFeatures)) {
        stop("Please define 'minFeatures' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$maxFeatures)) {
        stop("Please define 'maxFeatures' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$minCounts)) {
        stop("Please define 'minCounts' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$maxCounts)) {
        stop("Please define 'maxCounts' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$maxMito)) {
        stop("Please define 'maxMito' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      message("Running FilterDeadCells")
      seuratObject<-FilterDeadCells(seuratObject = seuratObject,
                                    minFeatures = optionalParameters$minFeatures,
                                    maxFeatures = optionalParameters$maxFeatures,
                                    minCounts = optionalParameters$minCounts,
                                    maxCounts = optionalParameters$maxCounts,
                                    maxMito = optionalParameters$maxMito)
    }
    else {
      tmp<-FilterDeadCellsByQuantile(seuratObject = seuratObject, lowQuantile = 0.1, highQuantile = 0.95)
      seuratObject<-tmp[[1]]
      svg(paste0(pathway,"QC_QuantileFiltering.svg"))
      print(tmp[[2]])
      dev.off()
      gg_preFiltering<-tmp[[2]]

    }

  }
  if (performNormalisation==TRUE) {
    seuratObject<-Seurat::NormalizeData(object = seuratObject,verbose = FALSE)
  }
  if(performVariableGeneDetection==TRUE){
    seuratObject<-Seurat::FindVariableFeatures(object = seuratObject, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  if (performScaling==TRUE) {
    seuratObject<-Seurat::ScaleData(object = seuratObject)
  }
  message("Imported ", length(seuratObject@meta.data$orig.ident), " cells from ", pathway, "with ID ", id, "\n")


  return(list(seuratObject, gg_preFiltering))
}



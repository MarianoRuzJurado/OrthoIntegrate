#' Summarize mapping statistics
#' This function reads in the summary file about the mapping with CellRanger or Starsolo and gives warnings if values surpass a threshold
#' @author Mariano Ruz Jurado
#' @param pathway A vector of pathways to the cellranger/star solo count result
#' @param id Vector of strings that are assigned to the concordant cells
#' @return a list with summary about mapping results
#' @export
SummarizeMapping <- function(pathway,id){
  ## mapping stats, i hope everyone uses cellranger and starsolo directories for these analysis, else no summary
  if (file.exists(paste0(unlist(strsplit(pathway,split = "outs"))[1],"outs/metrics_summary.csv"))) {
    metrics_summ.path <- paste0(unlist(strsplit(pathway,split = "outs"))[1],"outs/metrics_summary.csv")
    #define your own numeric class
    setClass('myNum')
    #define conversion
    setAs("character", "myNum",
          function(from) as.numeric(gsub(",","", gsub("%","",from))))
    #read data with custom colClasses
    metrics_summ <- read.csv(metrics_summ.path,
                             header = T,
                             stringsAsFactors=FALSE,
                             colClasses=c("myNum"))

    typeof(metrics_summ$Fraction.Reads.in.Cells)

    metrics_col <- as.data.frame(colnames(metrics_summ))
    rownames(metrics_col) <- metrics_col[,1]
    metrics_col[,1] <- as.character(as.vector(metrics_summ[1,]))
    metrics_summ <- metrics_col

    # warnings CELLRANGER
    if (metrics_summ[grep(pattern = "Confidently.to.Genome",rownames(metrics_summ)),] < 70) {
      warning(paste0(id,": Reads mapped confidently to genome only ", metrics_summ[grep(pattern = "Confidently.to.Genome",rownames(metrics_summ)),]))
    }
    if (metrics_summ[grep(pattern = "Confidently.to.Transcriptome",rownames(metrics_summ)),] < 30) {
      warning(paste0(id,": Reads mapped confidently to transcriptome only ", metrics_summ[grep(pattern = "Confidently.to.Transcriptome",rownames(metrics_summ)),]))
    }
    if (paste(unlist(strsplit(metrics_summ[grep(pattern = "Number.of.Cells",rownames(metrics_summ)),],split=",")),collapse = "") < 1000) {
      warning(paste0(id,": Estimated Number of Cells only ", metrics_summ[grep(pattern = "Number.of.Cells",rownames(metrics_summ)),]), " ,maybe the 0s were cut because of CR way of displaying numbers,  if unsure check CR web_summary")
    }
    if (as.numeric(paste(unlist(strsplit(metrics_summ[grep(pattern = "Median.Genes.per.Cell",rownames(metrics_summ)),],split=",")),collapse = "")) < 300) {
      warning(paste0(id,": Median Genes per Cell only ", metrics_summ[grep(pattern = "Median.Genes.per.Cell",rownames(metrics_summ)),])," ,maybe the 0s were cut because of CR way of displaying numbers, if unsure check CR web_summary")
    }
  } else {
    metrics_summ.path <- paste0(unlist(strsplit(pathway,split = "Gene"))[1],"Gene/Summary.csv")
    metrics_summ <- read.delim2(metrics_summ.path, header = F, sep = ",")
    rownames(metrics_summ) <- metrics_summ[,1]
    metrics_summ[,1] <- NULL

    # warnings STAR
    if (metrics_summ[7,] < 0.70) { # mapped to genome, no grep since same name as other row
      warning(paste0(id,": Reads mapped confidently to genome only ",metrics_summ[7,]))
    }
    if (metrics_summ[grep(pattern = "Transcriptome: Unique Genes",rownames(metrics_summ)),] < 0.30) {
      warning(paste0(id,": Reads mapped confidently to transcriptome only ", metrics_summ[grep(pattern = "Transcriptome: Unique Genes",rownames(metrics_summ)),]))
    }
    if (metrics_summ[grep(pattern = "Number of Cells",rownames(metrics_summ)),] < 1000) {
      warning(paste0(id,": Estimated Number of Cells only ", metrics_summ[grep(pattern = "Number of Cells",rownames(metrics_summ)),]))
    }
    if (metrics_summ[grep(pattern = "Median Genes per Cell",rownames(metrics_summ)),] < 300) {
      warning(paste0(id,": Median Genes per Cell only ", metrics_summ[grep(pattern = "Median Genes per Cell",rownames(metrics_summ)),]))
    }
  }
  return(metrics_summ)
}

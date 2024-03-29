% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FilterCells.R
\name{FilterDeadCells}
\alias{FilterDeadCells}
\title{Filter Seurat Object by setting values manually}
\usage{
FilterDeadCells(
  seuratObject,
  minFeatures = 300,
  maxFeatures = 6000,
  minCounts = 500,
  maxCounts = 15000,
  maxMito = 0.05
)
}
\arguments{
\item{seuratObject}{a Seurat object}

\item{minFeatures}{Number of minimal Features in a cell}

\item{maxFeatures}{Number of maximal Features in a cell}

\item{minCounts}{Number of minimal Counts in a cell}

\item{maxCounts}{Number of maximal Counts in a cell}

\item{maxMito}{Percentage of maximal mitochondrial content in a cell}
}
\value{
filtered seurat object
}
\description{
The number of features and UMIs (nFeature_RNA and nCount_RNA) are calculated for every object by Seurat.
For non-UMI data, nCount_RNA represents the sum of the non-normalized values within a cell
you can set the maximal percentage of mitochondrial features here and store it in object metadata as \code{percent.mito}.
We use raw count data since this represents non-transformed and non-log-normalized counts
The \% of UMI mapping to MT-features is a common scRNA-seq QC metric.
}
\author{
David John & Mariano Ruz Jurado
}

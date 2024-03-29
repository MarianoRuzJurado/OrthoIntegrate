% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OrthologueAssignment.R
\name{nucleotide.matching}
\alias{nucleotide.matching}
\title{nucleotide sequence comparison alignment}
\usage{
nucleotide.matching(
  mGene,
  replacement,
  OrthologueList_allHuman,
  species.1,
  species.2,
  alignment_type
)
}
\arguments{
\item{mGene}{Gene which will be matched against its possible orthologues in replacement}

\item{replacement}{vector containing possible orthologues found by biomaRt}

\item{OrthologueList_allHuman}{global orthologuelist which contains all Human genes and is filled with mouse orthologues}

\item{species.1}{name of the first species e.g. human, mouse, zebrafish, more to come...}

\item{species.2}{name of the second species e.g. human, mouse, zebrafish, more to come...}

\item{alignment_type}{Set the type of alignment for orthologue assignment, default is "global" for Needleman-Wunsch, use "local" for Smith-Waterman}
}
\description{
Reads in nucleotide sequence from NCBI for possible orthologues and calculates a score by pairwise alignment comparison
}
\author{
Mariano Ruz Jurado
}

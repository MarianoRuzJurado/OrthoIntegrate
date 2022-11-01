#' ensembl database orthologue settin using biomart function, converting human genes into mouse genes
#' @author Mariano Ruz Jurado
#' @param GeneNames Vector of Human gene names
#' @return List with converted Names based on ensembl database entries
#' @export
returnOrthologueslist <- function(GeneNames){
  human = biomaRt::useEnsembl("genes", host = "https://dec2021.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useEnsembl("genes", host = "https://dec2021.archive.ensembl.org", dataset = "mmusculus_gene_ensembl")
  #change attributes including gene_type, we need to know which homolog orthologues are real coding genes!!!
  genesV2 = biomaRt::getLDS(attributes = c("hgnc_symbol","description"), filters = c("hgnc_symbol"), values = GeneNames , mart = human, attributesL = c("mgi_symbol", "description"), martL = mouse, uniqueRows=T)
  return(genesV2)
}

#' function for uppercase the first letter of a character
#' @author Mariano Ruz Jurado
#' @param x character
#' @return character with first letter uppercase
#' @export
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#' Main function for Builduing a global useable Orthologue List
#'
#' This function uses GTF files for Human and Mice and starts building a Global useable Orthologue List
#' by assigning unique orthologues to human gene symbols. It uses Ensembl, Uniprot and NCBI databases.
#'
#'
#' @author Mariano Ruz Jurado
#' @param GTF.human path to GTF file for human genome
#' @param GTF.mice path to GTF file for mice genome
#' @return List with uniquely assigned Orthologues between Human and mice
#' @export
BuildOrthologues <- function(GTF.human, GTF.mice){

  #load in GTF human
  GTF.gr.h <- rtracklayer::import(GTF.human)
  GTF.df.h <- as.data.frame(GTF.gr.h)
  #load in GTF mice
  GTF.gr.m <- rtracklayer::import(GTF.mice)
  GTF.df.m <- as.data.frame(GTF.gr.m)

  #Build ensembl orthologuelist with human GTF
  httr::set_config(httr::config(ssl_verifypeer = F))
  OrthologueList_allHuman<-data.frame(HGNC.symbol=unique(GTF.df.h$gene_name), MouseGene=NA)
  df.orth<-returnOrthologueslist(unique(GTF.df.h$gene_name))

  #get all gtf mice names
  Mice.gtf.names<-data.frame(MGI.symbol=GTF.df.m$gene_name)

  #remove predicted genes from orthologue List
  df.orth <- df.orth[!grepl("predicted", df.orth$Gene.description.1),]
  df.orth <- df.orth[order(df.orth$HGNC.symbol),]

  #"SEPTIN" genes get a weird annotation of "Septin" in mouse genes, but they are annotated as "Sept" genes in the GTF mice -> change
  for (i in 1:length(df.orth$MGI.symbol)) {
    df.orth$MGI.symbol[i] <- gsub(pattern = "Septin", replacement = "Sept", df.orth$MGI.symbol[i])
  }

  #check the orthologues with the mice GTF, there are some strange ortholgue name founds which are not in the mice gtf
  df.orth <-df.orth[which(df.orth$MGI.symbol %in% Mice.gtf.names$MGI.symbol),] #lets kick them, we cant assign them correctly

  counter <- 0 # counts how often it fails to find an orthologue

  for (mGene in OrthologueList_allHuman$HGNC.symbol) {
    if(mGene %in% df.orth$HGNC.symbol){ #check if human gene has a mouse orthologue gene

      replacement <-df.orth[df.orth$HGNC.symbol==mGene,]$MGI.symbol
      if (length(replacement)==1 && replacement== '') { #check for empty ortholog
        replacement <- NA
      }
      if (length(replacement)==1 && replacement %in% OrthologueList_allHuman$MouseGene) { # dont let a mouse gene get another human ortholog, check if there is already the feature in the mouseGene list -> no double assignment
        replacement <- NA
      }

      if(length(replacement)>1){ #if more than one mouse gene represents the human gene

        #### SEQUENCE MATCHING Protein Sequence
        output.prot <- suppressWarnings(protein.matching(mGene,replacement,OrthologueList_allHuman)) ##list element 2 and 3 contains sequences
        OrthologueList_allHuman <- output.prot[[1]]

        ################ SEQUENCE MATCHING Nucleotide, if protein is empty
        if (is.null(output.prot[[2]]) || is.null(output.prot[[3]])) {

          output.nuc <- suppressWarnings(nucleotide.matching(mGene, replacement, OrthologueList_allHuman)) ##list element 2 contains sequence
          OrthologueList_allHuman <- output.nuc[[1]]

        }
        # if sequence matching doesnt give a hit because of no sequence entries in entrez(NCBI), make a match based on gene symbol similarity
        if (is.null(output.prot[[2]]) || is.null(output.prot[[3]]) && is.null(output.nuc[[2]])) {
          counter <- counter+1
          mGene.low <- tolower(mGene) # lowercase HGNC symbol for lower Levenshtein distances
          replacement <- replacement[replacement != ''] # delete empty character strings in vector, there is one feature
          replacement <- replacement[order(replacement)] # sort alphabetically
          Levensthein.DF <- RecordLinkage::levenshteinDist(mGene.low, replacement) # calculate Levenshtein distance for every hit
          Levensthein.DF <- sort(Levensthein.DF) # sort values increasing
          # , so in while loop there will be picked the right one
          low.score.hit <- which.min(Levensthein.DF) # take the lowest score by
          # taking values of column with lowest hit

          replacement.hit <- replacement[low.score.hit] # set it in variable
          j<-1 # set counter for while
          while (replacement.hit %in% OrthologueList_allHuman$MouseGene && !is.na(replacement.hit)) { # check if already in DF, if yes then take second hit
            replacement.hit <- replacement[low.score.hit+j]
            j<-j+1
          }
          OrthologueList_allHuman[OrthologueList_allHuman$HGNC.symbol==mGene,]$MouseGene=replacement.hit # set ortholog
        }
      }
      #only 1 replacement and not already in orthologue list, ideal match
      if (length(replacement)==1 && !(replacement %in% OrthologueList_allHuman$MouseGene))
      {
        OrthologueList_allHuman[OrthologueList_allHuman$HGNC.symbol==mGene,]$MouseGene=replacement
      }
    }
  }

  #get the uppercase overlaps into the orthologuelist
  MGI.symbols.up <- toupper(Mice.gtf.names$MGI.symbol)
  MGI.symbols.up <- unique(MGI.symbols.up)
  MGI.symbols.in <- MGI.symbols.up[MGI.symbols.up %in% OrthologueList_allHuman$HGNC.symbol == TRUE]
  #get rows of orthologue DT with missing mouse orthologue
  OrthologueList_allHuman.na <- OrthologueList_allHuman[rowSums(is.na(OrthologueList_allHuman))>0,]
  #check if the overlapping genes have matches with incomplete rows in the orthologue DF
  MGI.symbols.in.hits<- MGI.symbols.in[MGI.symbols.in %in% OrthologueList_allHuman.na$HGNC.symbol == TRUE]

  #set the found ortholgues in lower case into the orthologuelist
  for (ortholog.name in MGI.symbols.in.hits) {
    if (firstup(tolower(ortholog.name)) %in% OrthologueList_allHuman$MouseGene == FALSE){ # check if already in mouse Gene
      OrthologueList_allHuman[OrthologueList_allHuman$HGNC.symbol==ortholog.name,]$MouseGene=firstup(tolower(ortholog.name))
    }
  }

  #get the non NA rows
  OrthologueList <- OrthologueList_allHuman[complete.cases(OrthologueList_allHuman),]
  return(OrthologueList)
}

#' protein sequence matching alignment
#'
#' Reads in protein sequences from Uniprot for possible orthologues and calculates a score by pairwise alignment comparison
#'
#' @author Mariano Ruz Jurado
#' @param mGene  Gene which will be matched against its possible orthologues in replacement
#' @param replacement  vector containing possible orthologues found by biomaRt
#' @param OrthologueList_allHuman  global orthologuelist which contains all Human genes and is filled with mouse orthologues
#' @export
function (mGene, replacement, OrthologueList_allHuman)
{
  outh <- mygene::queryMany(mGene, scopes = "symbol", fields = c("entrezgene",
                                                                 "uniprot"), species = "human")
  outm <- try(mygene::queryMany(replacement, scopes = "symbol",
                                fields = c("entrezgene", "uniprot"), species = "mouse"),
              silent = T)
  if (!("try-error" %in% class(outm))) {
    uniprot.mice <- list()
    for (i in 1:nrow(outm)) {
      if (!is.null(unlist(outm[i, ]$uniprot.TrEMBL))) {
        uniprot.mice[i] <- outm[i, ]$uniprot.TrEMBL
      }
      else {
        uniprot.mice[i] <- NA
      }
      if (is.na(uniprot.mice[i]) && !is.null(unlist(outm[i, ]$uniprot.Swiss.Prot))) {
        uniprot.mice[i] <- outm[i, ]$uniprot.Swiss.Prot
      }
    }
    for (i in 1:length(uniprot.mice)) {
      if (!is.null(uniprot.mice[[i]][1]) && !is.na(uniprot.mice[[i]][1])) {
        uniprot.mice[[i]] <- uniprot.mice[[i]][1]
      }
      else {
        uniprot.mice[[i]] <- NA
        replacement[i] <- NA
      }
    }
    replacement <- replacement[!is.na(replacement)]
    uniprot.mice <- uniprot.mice[!is.na(uniprot.mice)]

    uniprot.human <- list()
    for (i in 1:nrow(outh)) {
      if (!is.null(unlist(outh[i, ]$uniprot.TrEMBL))) {
        uniprot.human[i] <- outh[i, ]$uniprot.TrEMBL
      }
      else {
        uniprot.human[i] <- NA
      }
      if (is.na(uniprot.human[i]) && !is.null(unlist(outh[i, ]$uniprot.Swiss.Prot))) {
        uniprot.human[i] <- outh[i, ]$uniprot.Swiss.Prot
      }
    }
    for (i in 1:length(uniprot.human)) {
      if (!is.null(uniprot.human[[i]][1]) && !is.na(uniprot.human[[i]][1])) {
        uniprot.human[[i]] <- uniprot.human[[i]][1]
      }
      else {
        uniprot.human[[i]] <- NA
      }
    }
    uniprot.human <- uniprot.human[!is.na(uniprot.human)]

    sequences.h <- suppressWarnings(UniprotR::GetSequences(unlist(uniprot.human),
                                                           directorypath = NULL))
    human.sequences <- sequences.h$Sequence

    sequences.m <- suppressWarnings(UniprotR::GetSequences(unlist(uniprot.mice),
                                                           directorypath = NULL))
    mice.sequences <- sequences.m$Sequence
    if (!is.null(human.sequences) && !is.null(mice.sequences)) {
      local.Align.list <- list()
      for (k in 1:length(mice.sequences)) {
        localAlign <- Biostrings::pairwiseAlignment(human.sequences,
                                                    mice.sequences[k], substitutionMatrix = "BLOSUM50",
                                                    gapOpening = 0, gapExtension = 8)
        local.Align.list[[k]] <- localAlign@score
      }
      names(local.Align.list) <- replacement
      local.Align.list <- local.Align.list[order(-unlist(local.Align.list))]
      replacement.hit <- names(local.Align.list[1])
      j <- 1
      while (replacement.hit %in% OrthologueList_allHuman$MouseGene &&
             !is.na(replacement.hit)) {
        replacement.hit <- names(local.Align.list[j +
                                                    1])
        j <- j + 1
      }
      OrthologueList_allHuman[OrthologueList_allHuman$HGNC.symbol ==
                                mGene, ]$MouseGene = replacement.hit
    }
  }
  else {
    human.sequences <- NA
    mice.sequences <- NA
  }
  return(list(OrthologueList_allHuman, human.sequences, mice.sequences))
}


#' nucleotide sequence comparison alignment
#'
#' Reads in nucleotide sequence from NCBI for possible orthologues and calculates a score by pairwise alignment comparison
#'
#' @author Mariano Ruz Jurado
#' @param mGene  Gene which will be matched against its possible orthologues in replacement
#' @param replacement  vector containing possible orthologues found by biomaRt
#' @param OrthologueList_allHuman  global orthologuelist which contains all Human genes and is filled with mouse orthologues
#' @export
nucleotide.matching <- function(mGene,replacement,OrthologueList_allHuman){


  out <- mygene::queryMany(mGene, scopes = "symbol", fields= c("entrezgene", "uniprot"),species ="human")
  outm <- try(mygene::queryMany(replacement, scopes = "symbol", fields= c("entrezgene", "uniprot"),species ="mouse"),silent = T) # use try for some genes it gives errors

  if (!("try-error" %in% class(outm))) {
    #human
    linked_seq_ids <- rentrez::entrez_link(dbfrom = "gene", id=out$entrezgene, db="nuccore")
    linked_transcripts <- linked_seq_ids$links$gene_nuccore_refseqrna
    head(linked_transcripts)
    if (!is.null(linked_transcripts)) { #check if it hits a sequence human

      all_recs <- rentrez::entrez_fetch(db="nuccore", id=linked_transcripts, rettype = "fasta")
      sequences <- unlist(strsplit(all_recs, split=">"))
      sequences.orth <- sequences[grep("NM",sequences)] # NM is for non predicted mRNAs
      if (identical(sequences.orth, character(0)) == T) {
        sequences.orth <- sequences[grep("NR", sequences)] ## check if it is non coding RNA if mRNA is 0
      }
      if (length(sequences.orth) > 1) {
        sequences.orth <- sequences.orth[grep("transcript variant 1", sequences.orth)]
      }


      sequences.orth <- unlist(strsplit(sequences.orth, split="\n"))
      Non.variant.h <- sequences.orth[2:(length(sequences.orth)-1)]
      Non.variant.h <- paste0(Non.variant.h, collapse="")
    }

    #mouse
    if (!is.null(outm$entrezgene)) {
      linked_seq_ids <- rentrez::entrez_link(dbfrom = "gene", id=outm$entrezgene, db="nuccore")
      linked_transcripts.m <- linked_seq_ids$links$gene_nuccore_refseqrna
      head(linked_transcripts.m)
    } else {
      linked_transcripts.m <- NULL
    }
    if (!is.null(linked_transcripts.m)) { #check if it hits a sequence mice

      all_recs <- rentrez::entrez_fetch(db="nuccore", id=linked_transcripts.m, rettype = "fasta")
      sequences <- unlist(strsplit(all_recs, split=">"))
      sequences.orth <- sequences[grep("NM",sequences)] # NM is for non predicted mRNAs
      if (identical(sequences.orth, character(0)) == T) {
        sequences.orth <- sequences[grep("NR", sequences)] ## check if there is non coding RNA if mRNA is 0
      }
      if (identical(sequences.orth, character(0)) == T) {
        sequences.orth <- sequences[grep("XM", sequences)] ## check if there is a predicted coding mRNA
      }
      if (identical(sequences.orth, character(0)) == T) {
        sequences.orth <- sequences[grep("XR", sequences)] ## check if there is a predicted non coding mRNA
      }

      orthologues.sequences <- as.list(sequences.orth)
      for (i in 1:length(orthologues.sequences)) {
        orthologues.seq <- unlist(strsplit(orthologues.sequences[[i]], split="\n"))
        orthologues.seq.2 <- orthologues.seq[2:length(orthologues.seq)]
        orthologues.sequences[[i]] <- paste0(orthologues.seq.2, collapse="")
        names(orthologues.sequences)[i] <- orthologues.seq[1]
      }
    }
    if (!is.null(linked_transcripts.m)) { # no transcript sequences, no alignment
      #alignment
      local.Align.list <- list()
      mat <- Biostrings::nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA")
      for (k in 1:length(orthologues.sequences)) {
        localAlign <- Biostrings::pairwiseAlignment(
          AAString(Non.variant.h),
          AAString(orthologues.sequences[[k]]),
          type="local",
          substitutionMatrix=mat , gapOpening = 5, gapExtension = 2)
        local.Align.list[[k]] <- localAlign@score
        #get gene name from string, stands between "()"
        names(local.Align.list)[k] <- stringr::str_match(names(orthologues.sequences)[k],"[(](.*?)[)]")[,2]
      }
      local.Align.list <- local.Align.list[order(-unlist(local.Align.list))] #order by score
      replacement.hit <- names(local.Align.list[1]) # first entry now has the highest score
      j<-1 # set counter for while
      while (replacement.hit %in% OrthologueList_allHuman$MouseGene && !is.na(replacement.hit)) { # check if already in DF, if yes then take second hit
        replacement.hit <- names(local.Align.list[j+1])
        j<-j+1
      }
      OrthologueList_allHuman[OrthologueList_allHuman$HGNC.symbol==mGene,]$MouseGene=replacement.hit # set orthologue
    }
  }
  else
    {
      linked_transcripts <- NA
    }
  return(list(OrthologueList_allHuman, linked_transcripts))
}

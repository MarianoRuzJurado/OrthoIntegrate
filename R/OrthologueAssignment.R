#' ensembl database orthologue setting using biomart function, converting mouse genes into human genes
#' @author Mariano Ruz Jurado
#' @param GeneNames Vector of Human gene names
#' @param species.1 name of the first species e.g. human, mouse, zebrafish, more to come...
#' @param species.2 name of the second species e.g. human, mouse, zebrafish, more to come...
#' @return List with converted Names based on ensembl database entries
#' @export
returnOrthologueslist <- function(GeneNames, species.1, species.2){
  human = biomaRt::useEnsembl("genes", host = "https://dec2021.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useEnsembl("genes", host = "https://dec2021.archive.ensembl.org", dataset = "mmusculus_gene_ensembl")
  zebrafish = biomaRt::useEnsembl("genes", host = "https://dec2021.archive.ensembl.org", dataset = "drerio_gene_ensembl")

  if (species.1 == "human") {
    mart.1 <- human
    attributes.1 <- c("hgnc_symbol","description")
    filters.1 <- c("hgnc_symbol")
  }
  if (species.1 == "mouse") {
    mart.1 <- mouse
    attributes.1 <- c("mgi_symbol","description")
    filters.1 <- c("mgi_symbol")
  }
  if (species.1 == "zebrafish") {
    mart.1 <- zebrafish
    attributes.1 <- c("external_gene_name","description")
    filters.1 <- c("")
  }
  if (species.2 == "human") {
    mart.2 <- human
    attributes.2 <-  c("hgnc_symbol","description")
  }
  if (species.2 == "mouse") {
    mart.2 <- mouse
    attributes.2 <- c("mgi_symbol","description")
  }
  if (species.2 == "zebrafish") {
    mart.2 <- zebrafish
    attributes.2 <- c("external_gene_name","description")
  }
  if (!(species.1 %in% c("human","mouse", "zebrafish")) || !(species.2 %in% c("human","mouse", "zebrafish")))  {
    stop("species.1 or species.2 is not one of 'human', 'mouse' or 'zebrafish'")
  }

  #change attributes including gene_type, we need to know which homolog orthologues are real coding genes!!!
  genesV2 = biomaRt::getLDS(attributes = attributes.1, filters = filters.1, values = GeneNames , mart = mart.1, attributesL = attributes.2, martL = mart.2, uniqueRows=T)
  return(genesV2)
}

#' function for uppercasing the first letter of a character
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
#' @param GTF.1 path to GTF file for e.g. human genome
#' @param GTF.2 path to GTF file for e.g. mice genome
#' @param species.1 name of the first species e.g. human, mouse, zebrafish, more to come...
#' @param species.2 name of the second species e.g. human, mouse, zebrafish, more to come...
#' @param alignment_type Set the type of alignment for orthologue assignment, default is "global" for Needleman-Wunsch, use "local" for Smith-Waterman
#' @return List with uniquely assigned Orthologues between Human and mice
#' @export
BuildOrthologues <- function(GTF.1, GTF.2, species.1, species.2, alignment_type="global"){

  if(!(alignment_type %in% c("global","local"))){
    stop("alignment type not defined as 'global' or 'local', see documentation for help.")
  }

  #load in GTF human
  GTF.gr.h <- rtracklayer::import(GTF.1)
  GTF.df.h <- as.data.frame(GTF.gr.h)
  #load in GTF mice
  GTF.gr.m <- rtracklayer::import(GTF.2)
  GTF.df.m <- as.data.frame(GTF.gr.m)

  #Build ensembl orthologuelist with human GTF
  httr::set_config(httr::config(ssl_verifypeer = F))
  OrthologueList_allHuman<-data.frame(HGNC.symbol=unique(GTF.df.h$gene_name), MouseGene=NA)
  df.orth<-returnOrthologueslist(unique(GTF.df.h$gene_name), species.1, species.2)

  #get all gtf mice names
  Mice.gtf.names<-data.frame(MGI.symbol=GTF.df.m$gene_name)

  #remove predicted genes from orthologue List
  df.orth <- df.orth[!grepl("predicted", df.orth$Gene.description.1),]
  df.orth <- df.orth[order(df.orth$HGNC.symbol),]
  #set colnames to ensure downstream compatability between species
  colnames(df.orth) <- c("HGNC.symbol","Gene.description","MGI.symbol","Gene.description.1")
  #"SEPTIN" genes get a weird annotation of "Septin" in mouse genes, but they are annotated as "Sept" genes in the GTF mice -> change
  for (i in 1:length(df.orth$MGI.symbol)) {
    df.orth$MGI.symbol[i] <- gsub(pattern = "Septin", replacement = "Sept", df.orth$MGI.symbol[i])
  }

  #check the orthologues with the mice GTF, there are some strange ortholgue name founds which are not in the mice gtf
  df.orth <-df.orth[which(df.orth$MGI.symbol %in% Mice.gtf.names$MGI.symbol),] #lets kick them, we cant assign them correctly

  counter <- 0 # counts how often it fails to find an orthologue

  for (mGene in OrthologueList_allHuman$HGNC.symbol) { # OrthologueList_allHuman$HGNC.symbol
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
        output.prot <- suppressWarnings(protein.matching(mGene,replacement,OrthologueList_allHuman,species.1,species.2, alignment_type)) ##list element 2 and 3 contains sequences
        OrthologueList_allHuman <- output.prot[[1]]

        ################ SEQUENCE MATCHING Nucleotide, if protein is empty
        if (is.null(output.prot[[2]]) || is.null(output.prot[[3]])) {

          output.nuc <- suppressWarnings(nucleotide.matching(mGene, replacement, OrthologueList_allHuman,species.1,species.2, alignment_type)) ##list element 2 contains sequence
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
  if (species.1 == "human") {
    MGI.symbols.up <- toupper(Mice.gtf.names$MGI.symbol)
    MGI.symbols.up <- unique(MGI.symbols.up)
    MGI.symbols.in <- MGI.symbols.up[MGI.symbols.up %in% OrthologueList_allHuman$HGNC.symbol == TRUE]
    #get rows of orthologue DT with missing mouse orthologue
    OrthologueList_allHuman.na <- OrthologueList_allHuman[rowSums(is.na(OrthologueList_allHuman))>0,]
    #check if the overlapping genes have matches with incomplete rows in the orthologue DF
    MGI.symbols.in.hits<- MGI.symbols.in[MGI.symbols.in %in% OrthologueList_allHuman.na$HGNC.symbol == TRUE]

  }

  #get the uppercase overlaps into the orthologuelist
  if (species.1 == "mouse") {
    MGI.symbols.up <- firstup(tolower(Mice.gtf.names$MGI.symbol))
    MGI.symbols.up <- unique(MGI.symbols.up)
    MGI.symbols.in <- MGI.symbols.up[MGI.symbols.up %in% OrthologueList_allHuman$HGNC.symbol == TRUE]
    #get rows of orthologue DT with missing mouse orthologue
    OrthologueList_allHuman.na <- OrthologueList_allHuman[rowSums(is.na(OrthologueList_allHuman))>0,]
    #check if the overlapping genes have matches with incomplete rows in the orthologue DF
    MGI.symbols.in.hits<- MGI.symbols.in[MGI.symbols.in %in% OrthologueList_allHuman.na$HGNC.symbol == TRUE]

  }

  #get the uppercase overlaps into the orthologuelist
  if (species.1 == "zebrafish") {
    MGI.symbols.up <- tolower(Mice.gtf.names$MGI.symbol)
    MGI.symbols.up <- unique(MGI.symbols.up)
    MGI.symbols.in <- MGI.symbols.up[MGI.symbols.up %in% OrthologueList_allHuman$HGNC.symbol == TRUE]
    #get rows of orthologue DT with missing mouse orthologue
    OrthologueList_allHuman.na <- OrthologueList_allHuman[rowSums(is.na(OrthologueList_allHuman))>0,]
    #check if the overlapping genes have matches with incomplete rows in the orthologue DF
    MGI.symbols.in.hits<- MGI.symbols.in[MGI.symbols.in %in% OrthologueList_allHuman.na$HGNC.symbol == TRUE]

  }

  #set the found ortholgues in lower case into the orthologuelist
  for (ortholog.name in MGI.symbols.in.hits) {

    if (species.2 == "human") {
      if (toupper(ortholog.name) %in% OrthologueList_allHuman$MouseGene == FALSE){ # check if already in species.1
        OrthologueList_allHuman[OrthologueList_allHuman$HGNC.symbol==ortholog.name,]$MouseGene=toupper(ortholog.name)
      }
    }

    if (species.2 == "mouse") {
      if (firstup(tolower(ortholog.name)) %in% OrthologueList_allHuman$MouseGene == FALSE){
        OrthologueList_allHuman[OrthologueList_allHuman$HGNC.symbol==ortholog.name,]$MouseGene=firstup(tolower(ortholog.name))
      }
    }

    if (species.2 == "zebrafish") {
      if (tolower(ortholog.name) %in% OrthologueList_allHuman$MouseGene == FALSE){
        OrthologueList_allHuman[OrthologueList_allHuman$HGNC.symbol==ortholog.name,]$MouseGene=tolower(ortholog.name)
      }
    }
  }

  #get the non NA rows
  OrthologueList <- OrthologueList_allHuman[complete.cases(OrthologueList_allHuman),]
  colnames(OrthologueList) <- c(species.1, species.2)
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
#' @param species.1 name of the first species e.g. human, mouse, zebrafish, more to come...
#' @param species.2 name of the second species e.g. human, mouse, zebrafish, more to come...
#' @param alignment_type Set the type of alignment for orthologue assignment, default is "global" for Needleman-Wunsch, use "local" for Smith-Waterman
#' @export
protein.matching <- function (mGene, replacement, OrthologueList_allHuman,species.1,species.2, alignment_type)
{
  outh <- try(mygene::queryMany(mGene, scopes = "symbol", fields= c("entrezgene", "uniprot"),species =species.1),silent = T)
  outm <- try(mygene::queryMany(replacement, scopes = "symbol",fields = c("entrezgene", "uniprot"), species = species.2),silent = T) #some mouse genes cant be handled by mygene
  if (!("try-error" %in% class(outm)) && !("try-error" %in% class(outh))) {
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

    sequences.h <- suppressWarnings(UniprotR::GetSequences(unlist(uniprot.human),directorypath = NULL)) # suppress the UniprotR warnings if not found Sequence
    human.sequences <- sequences.h$Sequence

    sequences.m <- suppressWarnings(UniprotR::GetSequences(unlist(uniprot.mice),directorypath = NULL)) # suppress the UniprotR warnings if not found Sequence
    mice.sequences <- sequences.m$Sequence

    if (!is.null(human.sequences) && !is.null(mice.sequences)) {
      local.Align.list <- list()
      for (k in 1:length(mice.sequences)) {

        #replace analogous Aminoacid with their common AC to be able to align
        if ("U" %in% strsplit(mice.sequences, split = "")[[k]]) {
          mice.sequences[k] <- paste(gsub("U", "C",strsplit(mice.sequences, split = "")[[k]]), collapse = "")
        }
        if ("U" %in% strsplit(human.sequences, split = "")[[1]]) {
          human.sequences <- paste(gsub("U", "C",strsplit(human.sequences, split = "")[[1]]), collapse = "")
        }
        #Run algorithm based on alignment type
        if (alignment_type == "global") {
          localAlign <- Biostrings::pairwiseAlignment(human.sequences,
                                                      mice.sequences[k], substitutionMatrix = "BLOSUM50",
                                                      gapOpening = 0, gapExtension = 8)
          local.Align.list[[k]] <- localAlign@score
        }
        if (alignment_type == "local") {
          localAlign <- Biostrings::pairwiseAlignment(human.sequences,
                                                      mice.sequences[k], substitutionMatrix = "BLOSUM50",
                                                      gapOpening = 0, gapExtension = 8, type="local")
          local.Align.list[[k]] <- localAlign@score
        }

      }
      names(local.Align.list) <- replacement
      local.Align.list <- local.Align.list[order(-unlist(local.Align.list))]
      replacement.hit <- names(local.Align.list[1])
      j <- 1
      while (replacement.hit %in% OrthologueList_allHuman$MouseGene && !is.na(replacement.hit)) { # check if already in DF, if yes then take second hit
        replacement.hit <- names(local.Align.list[j + 1])
        j <- j + 1
      }
      OrthologueList_allHuman[OrthologueList_allHuman$HGNC.symbol == mGene, ]$MouseGene = replacement.hit # set orthologue
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
#' @param species.1 name of the first species e.g. human, mouse, zebrafish, more to come...
#' @param species.2 name of the second species e.g. human, mouse, zebrafish, more to come...
#' @param alignment_type Set the type of alignment for orthologue assignment, default is "global" for Needleman-Wunsch, use "local" for Smith-Waterman
#' @export
nucleotide.matching <- function(mGene,replacement,OrthologueList_allHuman,species.1,species.2, alignment_type){


  out <- try(mygene::queryMany(mGene, scopes = "symbol", fields= c("entrezgene", "uniprot"),species =species.1),silent = T)
  outm <- try(mygene::queryMany(replacement, scopes = "symbol", fields= c("entrezgene", "uniprot"),species =species.2),silent = T) # use try for some genes there is no entry

  if (!("try-error" %in% class(outm)) && !("try-error" %in% class(out))) {
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
      Sys.sleep(0.1) # database could refuse to many requests
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

    if (!is.null(linked_transcripts.m) && !is.null(linked_transcripts)) { # no transcript sequences, no alignment
      #alignment
      local.Align.list <- list()
      mat <- Biostrings::nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA")
      #Run algorithm based on alignment type
      if (alignment_type == "global") {
        for (k in 1:length(orthologues.sequences)) {
          localAlign <- Biostrings::pairwiseAlignment(
            Biostrings::AAString(Non.variant.h),
            Biostrings::AAString(orthologues.sequences[[k]]),
            type="global",
            substitutionMatrix=mat , gapOpening = 5, gapExtension = 2)
          local.Align.list[[k]] <- localAlign@score
          #get gene name from string, stands between "()"
          names(local.Align.list)[k] <- stringr::str_match(names(orthologues.sequences)[k],"[(](.*?)[)]")[,2]
        }
      }

      if (alignment_type == "local") {
        for (k in 1:length(orthologues.sequences)) {
          localAlign <- Biostrings::pairwiseAlignment(
            Biostrings::AAString(Non.variant.h),
            Biostrings::AAString(orthologues.sequences[[k]]),
            type="local",
            substitutionMatrix=mat , gapOpening = 5, gapExtension = 2)
          local.Align.list[[k]] <- localAlign@score
          #get gene name from string, stands between "()"
          names(local.Align.list)[k] <- stringr::str_match(names(orthologues.sequences)[k],"[(](.*?)[)]")[,2]
        }
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

#' Subsets objects by orthologues found in the global list and gives Human gene names for Mouse objects
#' @author Mariano Ruz Jurado
#' @param SeuratObjectList.species.1 list with objects of species 1 to subset by found orthologues
#' @param OrthologueList previous mad elist with 1 to 1 orthologue assignment
#' @param SeuratObjectList.species.2 list with seurat objects of species 2 to subset by orthologues
#' @param species.1 name of the species, ("human","mouse", "zebrafish")
#' @param species.2 name of the species, ("human","mouse", "zebrafish")
#' @return list containing sublists with subsetted objects
#' @export
SubsetObjects <- function(OrthologueList,SeuratObjectList.species.1,SeuratObjectList.species.2, species.1, species.2){
  #for loop containing sub setting mice by all found orthologues and converting mice names in human orthologues
  SeuratObject.mouse.combined.orthologs.list <- list()
  human.converted <- list()
  #Ensure the colnames of OthologueList are set correctly
  colnames(OrthologueList) <- c(species.1, species.2)
  for (i in 1:length(SeuratObjectList.species.2)) {

    mouseGenes<-rownames(SeuratObjectList.species.2[[i]]@assays$RNA)
    mouseGenes.overlap <- character()
    human.names <- character()
    for (j in 1:length(mouseGenes)) {
      mGene <- mouseGenes[j]
      #check if the feature name is in the global ortholog list made of the human/mice gtf and uppercase matching
      if (mGene %in% OrthologueList[[species.2]] == TRUE){
        mouseGenes.overlap[j] <- mGene
      }
    }

    #Feature names with human ortholog
    mouseGenes.overlap <- mouseGenes.overlap[complete.cases(mouseGenes.overlap)]
    SeuratObject.mouse.combined.orthologs.list[[i]]<-subset(SeuratObjectList.species.2[[i]]
                                                            , features = mouseGenes.overlap)

    orthologlist.overlap <- rownames(SeuratObject.mouse.combined.orthologs.list[[i]])
    #get a list with converted mice to human gene names
    # , now that there are only mice features which have human ortholog
    for (k in 1:length(orthologlist.overlap)) {
      mGene <- orthologlist.overlap[k]
      human.names[k] <- OrthologueList[OrthologueList[[species.2]] == mGene,][[species.1]]
    }

    #insert the converted human.names into a list
    human.converted[[i]] <- human.names
    #get stats for sample
    stats.converting.mouse <- data.frame(rownames.length.object.before=length(mouseGenes)
                                   ,rownames.length.of.object.subsetted=length(orthologlist.overlap)
                                   ,rownames.length.object.converted=length(human.converted[[i]])
                                   ,row.names = levels(SeuratObjectList.species.2[[i]]$orig.ident))
    # #write out
    # write.xlsx(stats.converting,file = paste0(outputFolder,"/Orthologuestats/"
    #                                           ,levels(SeuratObjectList.mice[[i]]$orig.ident)
    #                                           ,"conversion.xlsx"), col.names = T,row.names = T)
  }

  #for loop containing subsetting human by all found orthologues
  SeuratObject.human.combined.orthologs.list <- list()
  for (i in 1:length(SeuratObjectList.species.1)) {

    humanGenes<-rownames(SeuratObjectList.species.1[[i]]@assays$RNA)
    humanGenes.overlap <- character()
    humanGenes.no.ortholog <- character()
    for (j in 1:length(humanGenes)) {
      hGene <- humanGenes[j]
      #check if the feature name is in the global ortholog list made of the human/mice gtf and uppercase matching
      if (hGene %in% OrthologueList[[species.1]] == TRUE){
        humanGenes.overlap[j] <- hGene
      }
    }

    #Feature names with human ortholog
    humanGenes.overlap <- humanGenes.overlap[complete.cases(humanGenes.overlap)]
    #Feature names without human ortholog
    humanGenes.no.ortholog <- humanGenes.no.ortholog[complete.cases(humanGenes.no.ortholog)]
    SeuratObject.human.combined.orthologs.list[[i]]<-subset(SeuratObjectList.species.1[[i]]
                                                            , features = humanGenes.overlap)
    #get stats for sample
    stats.converting.human <- data.frame(rownames.length.object.before=length(humanGenes)
                                   ,rownames.length.of.object.subsetted=length(humanGenes.overlap)
                                   ,row.names = levels(SeuratObjectList.species.1[[i]]$orig.ident))
    # #write out
    # write.xlsx(stats.converting,file = paste0(outputFolder,"/Orthologuestats/"
    #                                           ,levels(SeuratObjectList.human[[i]]$orig.ident)
    #                                           ,"conversion.xlsx"), col.names = T,row.names = T)

  }
  Resultlist <- list() # build a list with names
  Resultlist[["SeuratObject.species2.list"]] <- SeuratObject.mouse.combined.orthologs.list
  Resultlist[["SeuratObject.species1.list"]] <- SeuratObject.human.combined.orthologs.list
  Resultlist[["species1.converted.species2.names"]] <- human.converted
  return(Resultlist)
}

#' Renaming of Mice Objects with the found human names from SubsetObjects
#' @author Mariano Ruz Jurado
#' @param ObjList List with Seurat object to subset with orthologues names
#' @param newnames List with vector with new names generated by SubsetObjects function
#' @return Seurat object containing subsetted object with human orthologue names
#' @export
RenameGenesSeurat <- function(ObjList, newnames) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  SeuratObject.mouse.combined.orthologs.humanized.list <- list() #build return list
  for (i in 1:length(ObjList)) {
    obj <- ObjList[[i]]
    convert.names <- newnames[[i]]
    # print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
    RNA <- obj@assays$RNA

    if (nrow(RNA) == length(convert.names)) {
      if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- convert.names
      if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- convert.names
      # if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
    } else {"Unequal gene sets: nrow(RNA) != nrow(convert.names)"}
    obj@assays$RNA <- RNA
    SeuratObject.mouse.combined.orthologs.humanized.list[[i]] <- obj
  }
  return(SeuratObject.mouse.combined.orthologs.humanized.list)
}

#' Renaming of Objects with the found names from BuildOrthologues, works for Seuv5 objects
#' @author Mariano Ruz Jurado
#' @param ObjList List with Seuratv5 objects to subset with orthologues names
#' @param newnames List with vector with new names generated by SubsetObjects function
#' @return Seurat object containing subsetted object with human orthologue names, ready for Integration
#' @export
RenameGenesSeuratv5 <- function(ObjList, newnames) { # Replace gene names in different slots of a Seuratv5 object.
  SeuratObject.mouse.combined.orthologs.humanized.list <- list() #build return list
  for (i in 1:length(ObjList)) {
    obj <- ObjList[[i]]
    convert.names <- newnames[[i]]
    exp_mtx_counts <- as.matrix(obj@assays$RNA@counts)
    exp_mtx_data <- as.matrix(obj@assays$RNA@data)
    #data frame conversion
    con_df <- data.frame(species.1 = rownames(exp_mtx_counts),
                         species.2 = convert.names,
                         stringsAsFactors = F)
    # remove non fits, but shouldnt remove anything since its our orthoIntegrate 1:1 list
    con_df <- con_df[!is.na(con_df$species.2),,F]

    #subset by found genes
    exp_mtx_counts <- exp_mtx_counts[con_df$species.1,]
    exp_mtx_data <- exp_mtx_data[con_df$species.1,]

    #echange gene names
    rownames(exp_mtx_counts) <- con_df$species.2
    rownames(exp_mtx_data) <- con_df$species.2
    #Create new object with subsetted matrix
    obj_conv <- Seurat::CreateSeuratObject(counts = exp_mtx_counts, data = exp_mtx_data, meta.data = obj@meta.data )
    #Run find variable features for integration
    obj_conv <- Seurat::FindVariableFeatures(obj_conv, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    SeuratObject.mouse.combined.orthologs.humanized.list[[i]] <- obj_conv
  }
  return(SeuratObject.mouse.combined.orthologs.humanized.list)
}

#' Wrapper function for SubsetObject and RenameGenesSeurat as well as the Seurat integration steps
#' @author Mariano Ruz Jurado
#' @param SeuratObjectList.species.1 list with objects to subset by found orthologues
#' @param OrthologueList previously made list with 1 to 1 orthologue assignment
#' @param SeuratObjectList.species.2 list with seurat objects to subset by orthologues
#' @param species.1 name of the species, ("human","mouse", "zebrafish")
#' @param species.2 name of the species, ("human","mouse", "zebrafish")
#' @return Integrated Human/Mouse Seurat object with Human Nomenclature
#' @export
IntegrateObjects <- function(SeuratObjectList.species.1,SeuratObjectList.species.2,OrthologueList, species.1, species.2){
  SubsetList <- SubsetObjects(SeuratObjectList.species.1 = SeuratObjectList.species.1,
                              SeuratObjectList.species.2 = SeuratObjectList.species.2,
                              OrthologueList = OrthologueList,
                              species.1,
                              species.2)

  HumanizedList.mice <- RenameGenesSeurat(ObjList = SubsetList$SeuratObject.species2.list,
                                         newnames = SubsetList$species1.converted.species2.names)

  SeuratObjectList <- do.call("c",list(HumanizedList.mice,SubsetList$SeuratObject.species1.list))
  SeuratObject.anchors <- Seurat::FindIntegrationAnchors(object.list = SeuratObjectList, dims = 1:20)
  SeuratObject.combined <- Seurat::IntegrateData(anchorset = SeuratObject.anchors, dims = 1:20)

  return(SeuratObject.combined)
}

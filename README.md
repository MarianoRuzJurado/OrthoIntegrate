# OrthoIntegrate


<a href="https://zenodo.org/doi/10.5281/zenodo.10686730"><img src="https://zenodo.org/badge/559908015.svg" alt="DOI"></a>


R package functions for Orthologue Assignment and Integration of Single Cell Data between species

# <b> How to install </b>

Install through github:

```ruby
devtools::install_github("MarianoRuzJurado/OrthoIntegrate")
```

# <b> How to use it </b>

<b> Tested with Seuratv4, if you are using Seuratv5 please use the adjusted workflow as described below: </b>

First load in the ```OrthoIntegrate``` package:

```ruby
library(OrthoIntegrate)
```

Then we have to build a vector containing PATHS for our count matrices (e.g. 10X results: Directory containing barcodes.tsv, genes.tsv, matrix.mtx)
The Importer function needs also a vector of IDs for our samples:

```ruby
#example for human and mice data
resultList.human <- Importer(pathways = Sample.Paths.human,ids = Samplenames.human,
                             FilterCells = TRUE,
                             FilterByAbsoluteValues = FALSE,
                             performScaling = TRUE)
                             
resultList.mice <- Importer(pathways = Sample.Paths.mice,
                            ids = Samplenames.mice,
                            FilterCells = TRUE,
                            FilterByAbsoluteValues = FALSE,
                            performScaling = TRUE)
```
The function returns a list with the newly created seurat objects and a list of plots with information about the filtering.
If you dont want to use an automatically calculated filtering, you can set ```FilterByAbsoluteValues = TRUE``` and set thresholds manually:

```ruby
resultList.human <- Importer(pathways = Sample.Paths.human,
                             ids = Samplenames.human,
                             FilterCells = TRUE,
                             FilterByAbsoluteValues = TRUE,
                             performScaling = TRUE,
                             minFeatures=300,
                             maxFeatures=6000,
                             minCounts=500,
                             maxCounts=15000,
                             maxMito=0.05)
                             
resultList.mice <- Importer(pathways = Sample.Paths.mice,
                            ids = Samplenames.mice,
                            FilterCells = TRUE,
                            FilterByAbsoluteValues = TRUE,
                            performScaling = TRUE,
                            minFeatures=300,
                            maxFeatures=6000,
                            minCounts=500,
                            maxCounts=15000,
                            maxMito=0.05)
```

Additionally we can check our mapping statistics for the provided samples:

```ruby
# optional for summary of mapping results by CR or star solo
MappingSummary <- SummarizeMapping(pathways = Sample.Paths.human,
                                   ids = Samplenames.human) 
```

After our data is converted to seurat objects we may start to build a table containing orthologues for them. We will need this table to integrate samples from different species. The function needs GTF-files for our species:

```ruby
Orthologue.DF <- BuildOrthologues(GTF.1 = ".../Humangenes.gtf",
                                  GTF.2 = ".../Micegenes.gtf",
                                  species.1 = "human",
                                  species.2 = "mouse")

# I strongly recommend to save the Orthologue.DF file for future use.
```

It will start define orthologues for our genes by using the Ensembl, Uniprot and NCBI database by creating a global useable table.
After this step finished, we can subset our seurat objects by the found orthologues and integrate them into one object, <b> update for Seuratv5: </b> Please continue with ```SubsetObjects``` and ```RenameGenesSeuratv5``` function, since ```IntegrateObjects``` will not work with new layer structure:

```ruby
SeuratObject.combined <- IntegrateObjects(OrthologueList = Orthologue.DF,
                                          SeuratObjectList.species.1 = resultList.human$SeuratObjects,
                                          SeuratObjectList.species.2 = resultList.mice$SeuratObjects,
                                          species.1 = "human",
                                          species.2 = "mouse")
```
After succesful integration, we can continue our downstream analysis with an object containing single cell data from different species.

<hr>

If you are interested in how the objects changed after subsetting with the found orthologues use the ```SubsetObjects``` function:

```ruby
SubsetList <- SubsetObjects(OrthologueList = Orthologue.DF,
                            SeuratObjectList.species.1 = resultList.human$SeuratObjects,
                            SeuratObjectList.species.2 = resultList.mice$SeuratObjects,
                            species.1 = "human",
                            species.2 = "mouse")
```
This returns a list with the subsetted objects for our two species. Maybe you want to check the new nomeneclature on the objects. This can be achieved by using the ```RenamesGenesSeurat``` function on ```SubsetList```, <b> update for Seuratv5: </b> Please use ```RenamesGenesSeuratv5``` instead:

```ruby
HumanizedList.mice <- RenameGenesSeurat(ObjList = SubsetList$SeuratObject.species2.list,
                                        newnames = SubsetList$species1.converted.species2.names)
HumanizedList.mice <- RenameGenesSeuratv5(ObjList = SubsetList$SeuratObject.species2.list,
                                        newnames = SubsetList$species1.converted.species2.names)
```

if you used ```RenamesGenesSeuratv5``` then continue with:

```ruby
SeuratObjectList <- do.call("c",list(HumanizedList.mice, SubsetList$SeuratObject.species1.list))
SeuratObject.anchors <- FindIntegrationAnchors(object.list = SeuratObjectList)
SeuratObject.combined <- IntegrateData(anchorset = SeuratObject.anchors)

#You might need to rejoin the layers in the RNA assay before proceeding with further downstream analysis:
SeuratObject.combined <- JoinLayers(SeuratObject.combined, assay = "RNA")

```

# <b> Contribution Guidelines </b>
Raising up an issue in this Github repository might be the fastest way of submitting suggestions and bugs.
Alternatively you can write to my email: ruzjurado@med.uni-frankfurt.de

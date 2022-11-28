suppressMessages(library(ArchR))
library(here)
library(tidyverse)
library(parallel)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = min(length(fragments_fn), round(parallel::detectCores()/4)))
GENOMEDIR='/home/bnphan/resources/genomes/rn7'
load(file.path(GENOMEDIR,'rn7_liftoff_mm10NcbiRefSeq_ArchR_annotations.rda'))

#############################################################################
### make that arrow file from bgzipped fragments file (tsv.gz or bed.gz) ####
ArrowFiles <- here('data/raw_data/arrow') %>%
  list.files(pattern = '.arrow$', full.names = T)

## drop sample 6 b/c of low QC
ArrowFiles = ArrowFiles[!grepl('6C', ArrowFiles)]

###############################
### make the ArchR Project ####
dir.create(here('data/tidy_data/ArchRProjects'), showWarnings = F, recursive = T)
proj = ArchRProject( 
  ArrowFiles = ArrowFiles, copyArrows = TRUE,
  geneAnnotation = geneAnnotation, #must be the custom rn7 version
  genomeAnnotation = genomeAnnotation, #must be the custom rn7 version
  outputDirectory = here("data/tidy_data/ArchRProjects/Rat_Transgen_NAc_scATAC"))

## filter doublets
proj = filterDoublets( proj, cutEnrich = .5, cutScore = -Inf, filterRatio = 1)

## increase unique fragments cutoff to 10^3.5, remove cluster of low QC cell 
idxSample <- BiocGenerics::which(proj$nFrags > 10^3.6)
cellsSample <- proj$cellNames[idxSample]
proj = subsetCells(ArchRProj = proj, cellNames = cellsSample)

## Re-create the metadata from sample names
pd = data.frame(cellNames = getCellNames(proj)) %>% 
  mutate(sampleName = ss(cellNames, '#'), 
         Sire = case_when(grepl('C', sampleName) ~ 'Cocaine', 
                          grepl('M', sampleName) ~ 'Methamphetamine', 
                          TRUE ~ 'Saline'),
         ) %>% dplyr::select(-c(sampleName)) %>% 
  column_to_rownames('cellNames')

for(col in names(pd)){
  proj = addCellColData(
    ArchRProj = proj,
    data = pd[proj$cellNames,col],
    name = col,
    cells = getCellNames(proj),
    force = TRUE
  )
}
table(proj$Sample, proj$Sire)

## save project
proj = saveArchRProject(proj)


#################################
# do iterative LSI clustering 
varFeat = 200
pd = getCellColData(proj)
dimRed = names(attributes(proj)$reducedDims)
embedNames = names(attributes(proj)$embeddings)

iterLSIName = paste0("IterativeLSI",varFeat,'_ATAC')
print(iterLSIName)
if (iterLSIName %ni% dimRed){
  pdf()
  proj <- addIterativeLSI( proj, useMatrix = "TileMatrix", 
                           name = iterLSIName,
                           LSIMethod = 2, 
                           iterations = 4, # increase this if noticing subtle batch effects
                           scaleTo = 20000, # median unique fragment per cell
                           selectionMethod = 'var',
                           clusterParams = list( # See Seurat::FindClusters
                             resolution = c(.2, .5, .7, 1), # lower this if noticing subtle batch effects
                             sampleCells = 10000,  n.start = 10), 
                           varFeatures = varFeat * 1000, # also can reduce this if noticing subtle batch effects
                           dimsToUse = 1:30, force = FALSE)
  dev.off()
  proj = saveArchRProject(ArchRProj = proj)}


UMAPName = paste0("UMAPI",varFeat,'_ATAC')
if (UMAPName %ni% embedNames){
  print(UMAPName)
  proj <- addUMAP(proj, reducedDims = iterLSIName, 
                  name = UMAPName, nNeighbors = 30, minDist = 0.5, 
                  metric = "cosine", force = F)}


# add clusters
ClustersName = paste0("ClustersI",varFeat,'_ATAC')
if (ClustersName %ni% names(pd)){
  print(ClustersName)
  proj <- addClusters(proj, reducedDims = iterLSIName, method = "Seurat", 
                      algorithm = 2,
                      filterBias = TRUE, name = ClustersName, resolution = 1, force = T)
}

# add Harmony batch correction
HarmonyName = paste0("HarmonyI",varFeat,'_ATAC')
if (HarmonyName %ni% dimRed ){
  print(HarmonyName)
  proj <- addHarmony(proj, reducedDims = iterLSIName, 
                     max.iter.harmony = 15, name = HarmonyName, 
                     groupBy = c('Sample','Sire'), force = T)}


# add umap
UMAPName2 = paste0("UMAPH",varFeat,'_ATAC')
if (UMAPName2 %ni% embedNames){
  print(UMAPName2)
  proj <- addUMAP(proj, reducedDims = HarmonyName, 
                  name = UMAPName2, nNeighbors = 30, minDist = 0.5, 
                  metric = "cosine", force = F)}


# add clusters
ClustersName2 = paste0("ClustersH",varFeat,'_ATAC')
if (ClustersName2 %ni% names(pd)){
  print(ClustersName2)
  proj <- addClusters(proj, reducedDims = HarmonyName, method = "Seurat", 
                      algorithm = 2, 
                      name = ClustersName2, resolution = 1, force = T)}

proj = addImputeWeights(proj, reducedDims = HarmonyName)
proj = saveArchRProject(ArchRProj = proj)


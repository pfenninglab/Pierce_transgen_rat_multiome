suppressMessages(library(ArchR))
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

########################################
## load genomeAnnotation, geneAnnotation
GENOMEDIR='/pylon5/ibz3aep/bnphan/resources/genomes/rheMac10'
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))

#################################################
## create Arrow file and compute duplicate scores
for(PROJDIR in c('ArchR_snATAC_Striatum_MSN','ArchR_snATAC_Striatum_Other')){
  proj = loadArchRProject(path = here('Macaque_Multiome_Striatum/data/tidy_data/ArchRProjects',PROJDIR))
  for(varFeat in c(60)){
    pd = getCellColData(proj)
    dimRed = names(attributes(proj)$reducedDims)
    embedNames = names(attributes(proj)$embeddings)
    
    # add iterative LSI
    iterLSIName = paste0("IterativeLSIX",varFeat)
    if (iterLSIName %ni% dimRed){
      pdf()
      proj <- addIterativeLSI( proj, useMatrix = "TileMatrix", 
      name = iterLSIName,
      LSIMethod = 2, 
      iterations = 4, # increase this if noticing subtle batch effects
      scaleTo = 15000, # median unique fragment per cell
      selectionMethod = 'var',
      clusterParams = list( # See Seurat::FindClusters
        resolution = c(.1, .1, .3, .3), # lower this if noticing subtle batch effects
        sampleCells = 10000,  n.start = 10), 
      varFeatures = varFeat * 1000, # also can reduce this if noticing subtle batch effects
      dimsToUse = 1:30, force = TRUE)
      dev.off()
      proj = saveArchRProject(ArchRProj = proj)}
    
    # add Harmony batch correction
    HarmonyName = paste0("HarmonyX",varFeat)
    if (HarmonyName %ni% dimRed ){
      print(HarmonyName)
      proj <- addHarmony( proj, reducedDims = iterLSIName, 
                      max.iter.harmony = 20, name = HarmonyName, 
                      groupBy = c("Animal", 'Region', 'Batch'), force = TRUE)}
    
    # add umap for harmony
    UMAPName2 = paste0("UMAPX",varFeat)
    if (UMAPName2 %ni% embedNames){
      print(UMAPName2)
      proj <- addUMAP(proj, reducedDims = HarmonyName,  name = UMAPName2, 
                      nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)}
    
    # add clusters for Harmony
    ClustersName2 = paste0("ClustersX",varFeat)
    if (ClustersName2 %ni% names(pd)){
      print(ClustersName2)
      proj <- addClusters( proj, reducedDims = HarmonyName, 
                         method = "Seurat", name = ClustersName2, 
                         filterBias = TRUE,resolution = .5, force = TRUE)}
    
    proj = saveArchRProject(proj)
  }
}


####################################################################
### 2) re-cluster the main ArchR project after dropping the UNK cells
PROJDIR='Macaque_Multiome_Striatum/data/tidy_data/ArchRProjects'
ARCHDIR=here(PROJDIR,'ArchR_snATAC_Striatum')
proj = loadArchRProject(ARCHDIR)

# add iterative LSI
varFeat = 200
iterLSIName = paste0("IterativeLSI",varFeat,'_ATAC')
print(iterLSIName)
proj <- addIterativeLSI( proj, useMatrix = "TileMatrix", 
                         name = iterLSIName,
                         LSIMethod = 2, 
                         iterations = 4, # increase this if noticing subtle batch effects
                         scaleTo = 6000, # median unique fragment per cell
                         selectionMethod = 'var',
                         clusterParams = list( # See Seurat::FindClusters
                           resolution = c(.1, .1, .3, .3), # lower this if noticing subtle batch effects
                           sampleCells = 10000,  n.start = 10), 
                         varFeatures = varFeat * 1000, # also can reduce this if noticing subtle batch effects
                         dimsToUse = 1:30, force = TRUE)

proj = saveArchRProject(ArchRProj = proj)

# add Harmony batch correction
HarmonyName = paste0("HarmonyI",varFeat,'_ATAC')
print(HarmonyName)
proj <- addHarmony(proj, reducedDims = iterLSIName, 
                   groupBy = c("Animal", 'Region', 'Batch'),
                   max.iter.harmony = 30, name = HarmonyName, force = TRUE)

proj = addImputeWeights(proj, reducedDims = HarmonyName)

# add umap
UMAPName2 = paste0("UMAPH",varFeat,'_ATAC')
print(UMAPName2)
proj <- addUMAP(proj, reducedDims = HarmonyName, 
                name = UMAPName2, nNeighbors = 30, minDist = 0.5, 
                metric = "cosine", force = TRUE)

proj = saveArchRProject(ArchRProj = proj)

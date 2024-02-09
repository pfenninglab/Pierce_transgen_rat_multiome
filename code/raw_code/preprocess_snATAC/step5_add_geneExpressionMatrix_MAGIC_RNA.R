suppressMessages(library(ArchR))
library(SingleCellExperiment)
library(SeuratDisk)
library(Seurat)
library(SummarizedExperiment)
library(here); library(tidyverse)
library(rtracklayer)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = 8)

## read in gene annotations
annot = import('/home/bnphan/resources/genomes/rn7/rn7_liftoff_mm10_RefSeq.gtf')
names(annot) = annot$gene_name

########################################################
## load in the filtered MAGIC imputed snRNA-seq dataset 
obj = here::here('data/tidy_data/Seurat_projects',
                 "Rat_transgen_multiomeRNA_refined_all_SeuratObj_N5.h5Seurat") %>% 
  LoadH5Seurat(assay = 'MAGIC_RNA')

sce = as.SingleCellExperiment(obj)
sce = sce[!duplicated(rownames(sce)),]

## adding to ArchR needs to be integer counts
counts(sce) = round(2^logcounts(sce))
logcounts(sce) = NULL

# add seqnames to chr
sce = sce[rownames(sce) %in% names(annot), ]
rowRanges(sce) = annot[match(rownames(sce), names(annot))]
table(seqnames(rowRanges(sce)))

## trim cell barcodes to be correct
colnames(sce) = ss(colnames(sce), '_', 2)

## combine MAGIC RNA with the ArchR multiomeATAC cells
PROJDIR=here::here("data/tidy_data/ArchRProjects/Rat_Transgen_NAc_scATAC")
PROJDIR2=here::here("data/tidy_data/ArchRProjects/Rat_Transgen_NAc_multiome")
if(! dir.exists(PROJDIR2)){
  proj = loadArchRProject(path = PROJDIR)
  proj = saveArchRProject( ArchRProj = proj, outputDirectory = PROJDIR2)
} else{
  proj = loadArchRProject(path = PROJDIR2)
  proj2 = PROJDIR=here::here("data/tidy_data/ArchRProjects/Rat_Transgen_NAc_scATAC_clusterRat") %>% 
    loadArchRProject()
}

# 11190 cells have both ATAC and RNA profiles
table( proj$cellNames %in% colnames(sce))
table( colnames(sce) %in% proj$cellNames )

## Add the raw MAGIC RNA profiles 
# Overlap Per Sample w/ scATAC : 1S=262,2M=2764,3C=3076,4S=2308,5M=2780
proj <- addGeneExpressionMatrix(input = proj, seRNA = sce, force = TRUE)
proj <- proj[!is.na(proj$Gex_nUMI)]
proj = saveArchRProject(ArchRProj = proj, dropCells = T)

## transfer the peak set from one project to the other
proj = addPeakSet(proj, peakSet = getPeakSet(proj2))
proj <- addPeakMatrix(proj)
proj = saveArchRProject(ArchRProj = proj, dropCells = T)

###############################################
## add iterative LSI for the gene expression
pd = getCellColData(proj)
dimRed = names(attributes(proj)$reducedDims)
embedNames = names(attributes(proj)$embeddings)

iterLSIName = paste0("IterativeLSI_RNA")
print(iterLSIName)
if (iterLSIName %ni% dimRed){
  pdf()
  proj <- addIterativeLSI(
    ArchRProj = proj,
    clusterParams = list( resolution = 0.2,
                          sampleCells = 10000, n.start = 10 ),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix",
    depthCol = "Gex_nUMI", force = TRUE, 
    varFeatures = 2500, firstSelection = "variable",
    binarize = FALSE, name = iterLSIName)
  dev.off()
  proj = saveArchRProject(ArchRProj = proj)
}

# add Harmony batch correction
HarmonyName = paste0("HarmonyI_RNA")
if (HarmonyName %ni% dimRed ){
  print(HarmonyName)
  proj <- addHarmony(proj, reducedDims = iterLSIName,
                     max.iter.harmony = 15, name = HarmonyName,
                     groupBy = c('Sample', 'Sire'), force = T)
  }

# add umap
UMAPName2 = paste0("UMAPH_RNA")
if (UMAPName2 %ni% embedNames){
  print(UMAPName2)
  proj <- addUMAP(proj, reducedDims = HarmonyName,
                  name = UMAPName2, nNeighbors = 30, minDist = 0.5,
                  metric = "cosine", force = T)
  }

# add clusters
ClustersName2 = paste0("ClustersH_RNA")
if (ClustersName2 %ni% names(pd)){
  print(ClustersName2)
  proj <- addClusters(proj, reducedDims = HarmonyName, method = "Seurat",
                      name = ClustersName2, resolution = 1, force = T)
}
proj = saveArchRProject(ArchRProj = proj)


####################################################
## add combined ATAC + RNA dimensionality reductions
proj <- addCombinedDims(proj, reducedDims = c("HarmonyI200_ATAC", iterLSIName), name = "LSI_Combined")

# add Harmony batch correction
proj <- addHarmony(proj, reducedDims = 'LSI_Combined',
                   max.iter.harmony = 15, name = 'HarmonyI_Combined',
                   groupBy = c('Sample', 'Sire'), force = T)

# add umap
proj <- addUMAP(proj, reducedDims = 'HarmonyI_Combined',
                name = 'UMAPH_Combined', nNeighbors = 30, minDist = 0.5,
                metric = "cosine", force = T)

# add clusters
proj <- addClusters(proj, reducedDims = 'HarmonyI_Combined', method = "Seurat",
                    name = "ClustersH_Combined", resolution = 2, force = T)

proj = addImputeWeights(proj, reducedDims = "HarmonyI_Combined")
proj = saveArchRProject(ArchRProj = proj,)

# add peak2gene links matrix
proj <- addPeak2GeneLinks( ArchRProj = proj,dimsToUse = 1:30,
                           reducedDims = "HarmonyI_Combined", 
                           useMatrix = 'GeneExpressionMatrix',
                           scaleDims = TRUE, corCutOff = 0.75, k = 100, 
                           knnIteration = 500, overlapCutoff = 0.8, 
                           maxDist = 5e+05, scaleTo = 10^4, log2Norm = TRUE)

proj = saveArchRProject(ArchRProj = proj)




suppressMessages(library(ArchR))
library(ArchRtoSignac)
library(SeuratDisk)
library(Seurat)
library(Signac)
library(here)
library(tidyverse)
library(rtracklayer)
library(future)
library(harmony)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

#########################################################
# 0) Seurat uses the future package for parallelization
plan("multicore", workers = 12)
options(future.rng.onMisuse = 'ignore')
options(future.globals.maxSize = 60 * 1024^3)

##############################
## 1) load in the Seurat object
obj = here::here('data/tidy_data/Seurat_projects',
                 "Rat_transgen_multiomeRNA_refined_all_SeuratObj_N5.h5Seurat") %>% 
  LoadH5Seurat()

obj = RenameCells(obj, new.names = ss(colnames(obj), '_', 2))
head(colnames(obj))

##############################
## 2) load in the ArchR project
addArchRThreads(threads = 8)
PROJDIR=here::here("data/tidy_data/ArchRProjects/Rat_Transgen_NAc_scATAC_clusterRat")
proj = loadArchRProject(path = PROJDIR)

# 10959 cells have both ATAC and RNA profiles
table( proj$cellNames %in% colnames(obj))
table( colnames(obj) %in% proj$cellNames )

## subset ArchR project to cells in the 
proj <- proj[proj$cellNames %in% colnames(obj) ]
pkm <- getPeakMatrix(proj) 
names(proj@peakSet) = rownames(pkm) ## important for RunChromVar step
seqinfo(proj@peakSet) = seqinfo(BSgenome.Rnorvegicus.UCSC.rn7)[seqlevels(proj@peakSet)]
  
###########################################################################
## 3)  import Ensembl database from mouse and lift over to rn7 coordinates
library(EnsDb.Mmusculus.v79)
library(BSgenome.Rnorvegicus.UCSC.rn7)

annotations <- import('/home/bnphan/resources/genomes/rn7/rn7_liftoff_mm10_RefSeq.gtf')
seqinfo(annotations) = seqinfo(BSgenome.Rnorvegicus.UCSC.rn7)[seqlevels(annotations)]
annotations <- keepStandardChromosomes(annotations, pruning.mode = "coarse")
mcols(annotations)[,c('score', 'phase', 'source')] = NULL
table(annotations$type)

## convert ArchR to Signac 
# the directory before "/outs/" for all samples, need trailing /
fragments_dir <- "data/raw_data/fragments/" 

seurat_atac <- ArchR2Signac(
  ArchRProject = proj,
  refversion = "rn7",
  fragments_dir = fragments_dir,
  pm = pkm, # peak matrix from getPeakMatrix()
  ## put in fragments_dir+ sample + "outs", "fragments.tsv.gz"
  fragments_fromcellranger = "Y", 
  fragments_file_extension = '.tsv.gz', 
  annotation = annotations # annotation from getAnnotation()
)

###############################################################################
## 4) stash the Signac data in the multiome RNA object and process ATAC data
head(colnames(seurat_atac))
head(colnames(obj))
obj = RenameCells(obj, new.names = gsub('#', '_', colnames(obj)))
obj = obj[,colnames(seurat_atac)] ## subset to just the RNA cells w/ ATAC data

## add the Signac Peak object to the Seurat object w/ multiomeRNA
obj[['ATAC']] = seurat_atac[['peaks']]

## add dimensionality reduction for peaks w/ Signac
DefaultAssay(obj) <- 'ATAC'
obj <- obj %>% FindTopFeatures(min.cutoff = 5) %>% 
  RunTFIDF() %>% RunSVD() %>% 
  RunHarmony(reduction = "lsi", group.by.vars = 'Sample',  dims.use = 2:40, 
             assay.use = "ATAC", project.dim = F, reduction.save = "harmony.lsi")

Idents(obj) <- "cluster_rat"
levels(obj) <- c("D1-MSN", "D1/D3-MSN", "D2-MSN", "D1/D2-Hybrid-MSN", 
                  "D1-NUDAP-MSN", "D1-ICj-MSN", "Pvalb-INT", "Sst-INT" ,
                  "Astrocytes", "Microglia", "Oligos", "Oligos_Pre")
pdf()
DepthCor(obj, reduction = "harmony.lsi")
dev.off()

## the WNN multimodal join NN building method
obj <- FindMultiModalNeighbors(
  obj, reduction.list = list("harmony", "harmony.lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight")

# build a joint UMAP visualization
obj <- RunUMAP( object = obj, nn.name = "weighted.nn", assay = "RNA")

## compute the GC bias background
obj <- RegionStats(obj, genome = BSgenome.Rnorvegicus.UCSC.rn7, assay = 'ATAC')

###########################################################
## 5) compute motif annotations and chromVar deviations
library(motifmatchr)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(TFutils)

## use human set b/c mouse set is incomplete, 
## https://github.com/stuart-lab/signac/issues/58#issuecomment-571183108
pfm <- getMatrixSet( x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information to peaks
obj <- AddMotifs(obj, genome = BSgenome.Rnorvegicus.UCSC.rn7, pfm = pfm)
obj <- RunChromVAR(object = obj, genome = BSgenome.Rnorvegicus.UCSC.rn7, assay = 'ATAC')

obj[['ATAC']]@motifs
names(attributes(obj[['ATAC']]@motifs))

## set this so H5Seurat writing doesn't crash
obj[['ATAC']]@motifs@meta.data = obj[['ATAC']]@motifs@motif.names %>% 
  unlist() %>% enframe() %>% dplyr::rename('motif' = 'name', 'motif.names' = value)

#####################################
## 6) save the Seurat/Signac object
SaveH5Seurat(obj, here::here('data/tidy_data/Seurat_projects',
           "Rat_transgen_RNA-ATAC_SeuratObj_N5.h5Seurat"),overwrite =TRUE)


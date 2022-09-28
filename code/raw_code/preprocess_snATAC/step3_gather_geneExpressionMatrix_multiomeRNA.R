library(scds) ## doublet removal
library(SingleCellExperiment) ## works w/ SCDS
library(SeuratDisk)
library(Seurat)
library(Matrix)
library(here)
library(tidyverse)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

# read in srna-seq matrix files, these folders should have a barcodes.tsv, genes.tsv, and matrix.mtx files
STARsoloDIR= here('data/raw_data/CellRanger_outs')
STARsolo_fn = file.path(list.dirs(STARsoloDIR, recursive = FALSE), 'filtered_feature_bc_matrix')
names(STARsolo_fn) = ss(STARsolo_fn, '/', 9)
STARsolo_fn

## read in GEM wells and compute scds doublet scores
countList <- STARsolo_fn %>% lapply(Read10X)
countList <- lapply(countList, function(x) x[["Gene Expression"]])
seuratList = mapply(CreateSeuratObject, counts = countList, project =names(countList))
sceList = lapply(seuratList, as.SingleCellExperiment)

## filter out inferred doublets
sceList <- sceList %>% lapply(cxds_bcds_hybrid, list("retRes"=TRUE))
sceList <- lapply(sceList, cxds_bcds_hybrid,list("retRes"=TRUE))
sceList <- lapply(sceList, function(sce) {
  mcols(sce) = mcols(sce)[,!names(mcols(sce)) %in% c('cxds_hvg_bool', 'cxds_hvg_ordr')]
  sce = sce[, colData(sce)$hybrid_score < 1.0 ]
  return(sce)
})

## combine GEM wells
sce = do.call('cbind', sceList)
rownames(sce) = rownames(sceList[[1]])
sce$Sample = sce$orig.ident
sce$Barcode = paste0(sce$Sample, '#', colnames(sce))
sce$Sire = ifelse(grepl('M', sce$Sample), 'Methamphetamine', 
                  ifelse(grepl('C', sce$Sample), 'Cocaine', 'Saline'))
colnames(sce) = sce$Barcode

## save 1 object w/ the raw counts from the multiomeRNA side of things
sce_fn = file.path('data/tidy_data/rdas', 'raw_multiomeRNA_transgen_rat.sce.rds')
saveRDS(sce, sce_fn)

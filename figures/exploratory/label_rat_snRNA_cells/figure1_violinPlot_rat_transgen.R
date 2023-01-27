## packages for data table processing 
library(here)
library(tidyverse)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

FIGDIR='figures/exploratory/label_rat_snRNA_cells'
dir.create(here(FIGDIR, 'plots'), recursive = T, showWarnings = F)

######################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 12 cores
plan("multicore", workers = 4)
options(future.globals.maxSize = 100 * 1024^3)
options(future.rng.onMisuse = 'ignore')

celltypes = c('D1-MSN','D1/D3-MSN','D2-MSN','D1-ICj-MSN','D1-NUDAP-MSN','D1/D2-Hybrid-MSN', 
              'Pvalb-INT','Sst-INT', 'Astrocytes', 'Microglia', 'Oligos', 'Oligos_Pre')

######################################################
# 1) Seurat uses the future package for parallelization
## load in the full dataset for label refinement
obj = here('data/tidy_data/Seurat_projects', 
           "Rat_transgen_multiomeRNA_refined_all_SeuratObj_N5.h5Seurat") %>% 
  LoadH5Seurat()
Idents(obj) = 'cluster_rat'
table(Idents(obj) %in%celltypes )
levels(obj) = celltypes

## load in the MSN dataset for label refinement
obj_msn = here('data/tidy_data/Seurat_projects', 
                       "Rat_transgen_multiomeRNA_refined_msn_SeuratObj_N5.h5Seurat") %>% 
  LoadH5Seurat()

cols_rat = ArchR::paletteDiscrete(unique(obj$cluster_rat))
cols_mac = ArchR::paletteDiscrete(unique(obj$cluster_macaque))

######################################################
# 1) Seurat uses the future package for parallelization
width = 8; height = 4; 

fn = here(FIGDIR, 'plots', 'rat_transgen_snRNA_celltype_markers.all.png')
png(fn, units = 'in',  res = 300, width = width*2, height = height*2)
features =  c('Rbfox3', 'Lhx6', 'Cx3cr1', 'Mog', 'Pdgfra', 'Aqp4')
VlnPlot(obj, features = features,  group.by = 'cluster_rat',  cols = cols_rat,
        assay = 'MAGIC_RNA', pt.size = 0) & 
  theme(axis.title.x=element_blank())
dev.off()


fn = here(FIGDIR, 'plots', 'rat_transgen_snRNA_interneuron_markers.all.png')
png(fn, units = 'in',  res = 300, width = width*2, height = height*2)
features2 = c('Pvalb', 'Pthlh', 'Chrna3','Sst', 'Npy',  'Chat')
VlnPlot(obj, features = features2,  group.by = 'cluster_rat', cols = cols_rat,
        assay = 'MAGIC_RNA', pt.size = 0) & 
  theme(axis.title.x=element_blank())
dev.off()


fn = here(FIGDIR, 'plots', 'rat_transgen_snRNA_msn_markers.all.png')
png(fn, units = 'in',  res = 300, width = 6*2, height = height*2)
features =  c('Drd1', 'Drd2', 'Drd3', 'Cpne4', 'Rxfp1', 'Oprm1')
VlnPlot(obj_msn, features = features, group.by = 'cluster_rat',  cols = cols_rat,
        assay = 'MAGIC_RNA', pt.size = 0) & 
  theme(axis.title.x=element_blank())
dev.off()

fn = here(FIGDIR, 'plots', 'rat_transgen_snRNA_msn_umap.all.png')
png(fn, units = 'in',  res = 300, width = 2*2, height = 3*2)
DimPlot(obj_msn, reduction = "umap", group.by = "cluster_rat", label = TRUE, label.size = 5, 
        cols = cols_rat) & theme(legend.position = 'none')
dev.off()


VlnPlot(obj_merged, features = features2, group.by = 'cluster_rat', assay = 'MAGIC_RNA', pt.size = 0)







DimPlot(obj_msn, reduction = "umap", group.by = "cluster_rat", label = TRUE, 
        label.size = 10, cols = cols_rat)

DimPlot(obj_msn, reduction = "umap", group.by = "cluster_rat", label = TRUE, 
        label.size = 10, cols = cols_rat)





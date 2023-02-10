library(readxl)
library(tidyverse)
library(here)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)
library(BiocParallel)
library(doParallel)

DATADIR = 'data/tidy_data/scenic_grn'
dir.create(here(DATADIR,'loom'), recursive = T, showWarnings = F)

## load in the filtered Seurat object
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "Rat_transgen_multiomeRNA_refined_all_SeuratObj_N5.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')

DefaultAssay(object = obj_merged) = 'RNA'
names(obj_merged[[]])
table(obj_merged$Sample)
table(obj_merged$cluster_rat) %>% sort()

## downsample oligos to 2nd most populous cell type, astrocytes
## this reduces compute by 13%
ind_drop = which(obj_merged$cluster_rat == 'Oligos') %>% sample(size = 2489)
obj_subset = obj_merged[, -ind_drop]

## subset genes to just those analyzed in differential analyses
## this reduces compute by 33%
res.celltype = here('data/tidy_data/differential_expression/rdas', 
                    'diff_gene_3-way_countsplit.rds') %>% readRDS()
genes = res.celltype %>% pull(gene) %>% sort() %>% unique()
obj_subset2 = obj_subset[genes, ]
table(obj_subset2$Sample)
table(obj_subset2$cluster_rat) %>% sort()

## create loom object for each sample
for( sub in obj_merged$Sample %>% sort() %>% unique()){
  out_loom_fn = here(DATADIR,'loom', paste0("Rat_transgen_multiomeRNA_refined.",sub,".loom"))
  if(!file.exists(out_loom_fn)){
    print(paste('making loom for:', sub))
    obj_subset3 = obj_subset2 %>% subset(Sample == sub)
    out_loom <- as.loom(obj_subset3, filename = out_loom_fn, verbose = T)
    out_loom$close_all() ## close 
  }
}


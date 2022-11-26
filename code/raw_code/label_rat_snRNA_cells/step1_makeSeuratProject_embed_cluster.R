## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(harmony)
library(SeuratDisk)
library(SeuratWrappers)
library(SingleCellExperiment)

library(future)
library(Rmagic)
library(phateR)
library(flexmix)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

###################################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 28 cores with 500Gb
plan("multicore", workers = 24)
options(future.rng.onMisuse = 'ignore')
options(future.globals.maxSize = 150 * 1024^3)

###########################################################################
# 1) load in indvidual snRNA-seq objects and create merged Seurat projects
obj = here('data/tidy_data/rdas/raw_multiomeRNA_transgen_rat.sce.rds') %>% 
  readRDS() %>% as.Seurat()

## add % ribosomal or mitochondrial genes
obj[["percent.ribo.mito"]] <- PercentageFeatureSet(object = obj, pattern = "^Rp[sl]|^Mt")
obj$scds.keep = ifelse(obj$hybrid_score < 1.0, 'keep', 'doublet')

## drop the low QC sample
obj = obj[, obj$Sample != '6C']
obj$Sample = droplevels(obj$Sample)
num_samples = length(unique(obj$Sample))

## compute miQC quality scores per sample
objList = SplitObject(obj, split.by = "Sample") %>% 
  ## estimate cells w/ high score likely "compromised"/low QC cells,
  lapply(RunMiQC, percent.mt = "percent.ribo.mito", nFeature_RNA = "nFeature_RNA", 
         posterior.cutoff = 0.75, model.slot = "flexmix_model") 

obj = merge(objList[[1]], objList[-1], names(objList))

###########################################################################
# 2) integrate the samples into each other and remove batch effects in UMAP
obj <- obj %>% 
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.ribo.mito") %>% 
  ## use Harmony to remove batch effects across samples
  RunPCA(verbose = TRUE) %>% 
  RunHarmony(group.by.vars = 'Sample', assay.use = "SCT", verbose = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  ## create data-based clustering
  FindNeighbors(dims = 1:30, verbose = TRUE) %>% 
  FindClusters(resolution = 1, algorithm = 2, verbose = TRUE)


###########################################################################
# 3) find the low QC clusters based on doublets and miQC estimates
## look at which clusters should be kept by miQC per-cell fraction
obj@meta.data %>% group_by(seurat_clusters) %>%
  summarise(num = n(), prop = sum(miQC.keep == 'keep') / n() ) %>% 
  arrange(prop)

## look at which clusters should be kept by both metrics
(t1 = obj@meta.data %>% group_by(seurat_clusters) %>%
    summarise(num = n(), 
              numKeep = sum(scds.keep == 'keep' & miQC.keep == 'keep'), 
              prop = sum(numKeep) / n() ) %>% arrange(prop))

## keep cells in the clusters that have more than 10% of OK cells
good_clusters <- t1 %>% filter(prop > 0.10) %>% pull(seurat_clusters)

## export unfiltered per-cell QC metric table
dir.create(here('data/tidy_data/tables'), showWarnings = F)
save_qcTale_fn = here('data/tidy_data/tables', 
                      paste0("Rat_transgen_multiomeRNA_unfiltered_QC_table_N",num_samples,'.txt.gz'))
write_tsv(obj@meta.data, save_qcTale_fn)

## subset cells to those not predicted low QC or doublet
obj_filtered = subset(obj, subset = miQC.keep == "keep" & scds.keep == "keep" & 
                        seurat_clusters %in% good_clusters)

## recompute PCA and UMAP embedding post-filtering
obj_filtered = obj_filtered %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>% 
  FindClusters(resolution = 0.5, algorithm = 2, verbose = FALSE)

#################################################################################
# 4) compute the MAGIC imputed matrices for both the RNA and SCT normalized objects
DefaultAssay(obj_filtered) <- "SCT"
obj_filtered <- magic(obj_filtered, genes='all_genes')

## smooth out the genes with imputation across full data
DefaultAssay(obj_filtered) <- "RNA"
obj_filtered <- obj_filtered %>% 
  NormalizeData() %>% magic(genes='all_genes')

## change back to the MAGIC_SCT imputed matrix
DefaultAssay(obj_filtered) <- "SCT"


#########################################################
# 6) save projects, as h5 object for on-disk computation
dir.create(here('data/tidy_data/Seurat_projects'), showWarnings = F)
save_proj_h5_fn = here('data/tidy_data/Seurat_projects', 
                       paste0("Rat_transgen_multiomeRNA_filteredImputed_SeuratObj_N",num_samples,'.h5Seurat'))
SaveH5Seurat(obj_filtered, filename = save_proj_h5_fn,overwrite = TRUE)


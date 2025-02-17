## packages for further processing 
library(here)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(DropletQC)
library(future)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

###################################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 28 cores with 500Gb
plan("multicore", workers = 1)
options(future.rng.onMisuse = 'ignore')
options(future.globals.maxSize = 150 * 1024^3)

###########################################################################
# 1) load in indvidual snRNA-seq objects and create merged Seurat projects
save_fn = list.files(here('/projects/pfenninggroup/singleCell/Jarvis_songbird_snRNA-seq/data/raw_data/Seurat_objects_combined'), 
                     pattern = 'STARsolo_SoupX_rawCounts', full.names = T)
names(save_fn) = basename(save_fn) %>% ss('STARsolo_SoupX_rawCounts_|.rds', 2)
num_samples = length(save_fn)
objList = lapply(save_fn, readRDS)

########################################################
# 2) use Seurat reciprocal PCA to join samples together
## find integrating features
features <- SelectIntegrationFeatures(object.list = objList, nfeatures = 3000)
objList <- PrepSCTIntegration(object.list = objList, anchor.features = features)
objList <- lapply(X = objList, FUN = RunPCA, features = features, verbose = FALSE)

## select representative caudate and putamen samples for integration
## both these samples are control subjects/good number cells & depth
## and one of each sex based on:
## https://satijalab.org/seurat/articles/integration_large_datasets.html
ref = which(names(objList) %in% c('AN1', 'AS1', 'ASH1', 'ASS1', 'AX1', 'AXH1', 'AXS1', 'DB1', 'DH1', 'DS1', 'HB1', 'HH1', 'HSS2', 'LAI2', 'LAIH1', 'LS1', 'LB', 'LH1', 'LMAN1', 'LS1', 'LSH1', 'LSS1', 'PH1', 'PS1', 'RAH1', 'RAS1'))

## find pair-wise anchoring cells between samples and each reference
obj.anchors <- FindIntegrationAnchors(
  object.list = objList, normalization.method = "SCT", reference = ref,
  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)

## merging samples together into a joint space
#obj_merged <- obj.anchors %>% 
obj_merged <- IntegrateData(obj.anchors, normalization.method = "SCT", dims = 1:30, k.weight = 90)
obj_merged <- RunPCA(obj_merged, verbose = FALSE) 
obj_merged <- RunUMAP(obj_merged, reduction = "pca", dims = 1:30)
obj_merged <-  FindNeighbors(obj_merged, ims = 1:30, verbose = TRUE)
obj_merged <-  FindClusters(obj_merged, resolution = 1, algorithm = 2, verbose = TRUE)


#####################################
# 3) add in patient/sample metadata
#pheno = readxl::read_xlsx("/projects/pfenninggroup/singleCell/Jarvis_songbird_snRNA-seq/data/tidy_data/tables/finch_snRNAseq_samples_metadata.xlsx")

#head(pheno)
# undo previous adds
#obj_merged@meta.data = obj_merged@meta.data[,1:13]
#ind = match(obj_merged$orig.ident, pheno$`Bird ID`)

# find the cells that don't have a meta data entry
#table(obj_merged$orig.ident[!obj_merged$orig.ident %in%  pheno$`Bird ID`])
#
# use indexing of phenotype data to add meta
#obj_merged@meta.data = cbind(obj_merged@meta.data, pheno[ind,])

#head(obj_merged@meta.data)

# check all cells have metadata column, should say TRUE
#all.equal(obj_merged$orig.ident, singing$`Bird ID`)

#pheno = readxl::read_xlsx(here('data/tidy_data/tables/finch_snRNAseq_samples_metadata.xlsx')) %>%
#  column_to_rownames('sampleid')

## take a look at the phenotype table
#head(pheno)

## look at the per-cell meta data
#head(obj_merged@meta.data)

## append pt. phenotypes to single cell meta data 
#obj_merged@meta.data = cbind(obj_merged@meta.data, pheno[obj_merged[[]][,'orig.ident'],])
#head(obj_merged@meta.data)

## this section was not run due to lack of genefull/velocity data

###############################################
# 5) filter the lower quality cells

## Run DropletQC's identify_empty_drops function
#nf.umi <- obj_merged[[]] %>% mutate(nf=dropletQC.nucFrac, umi=nCount_RNA) %>% 
#  relocate(nf, umi, .before = everything())

## estimate cell, empty droplet, and damaged cell w/ DropletQC algorithms
#DropletQC.ed <- nf.umi  %>% 
#  identify_empty_drops(nf_rescue = 0.50, umi_rescue = 1000) %>%
#  relocate(cell_status, seurat_clusters, .after = umi) %>%
#  dplyr::select(nf:seurat_clusters) %>%
#  identify_damaged_cells(nf_sep = 0.15, umi_sep_perc = 50, verbose = FALSE)
#DropletQC.ed = DropletQC.ed$df

## add Droplet QC empty droplet estimation to metadata
#obj_merged$dropletQC.keep = DropletQC.ed[colnames(obj_merged), 'cell_status']

## look at which clusters should be kept by miQC per-cell fraction
obj_merged@meta.data %>% group_by(seurat_clusters) %>%
  summarise(num = n(), prop = sum(miQC.keep == 'keep') / n() ) %>% 
  arrange(prop)

## look at which clusters should be kept by doublet SCDS per-cell fraction
obj_merged@meta.data %>% group_by(seurat_clusters) %>%
  summarise(num = n(), prop = sum(scds.keep == 'keep') / n() ) %>% 
  arrange(prop)

## look at which clusters should be kept by both metrics
(t1 = obj_merged@meta.data %>% group_by(seurat_clusters) %>%
    summarise(num = n(), 
              numKeep = sum(scds.keep == 'keep' & miQC.keep == 'keep'), # & dropletQC.keep == 'cell'), 
              prop = sum(numKeep) / n() ) %>% arrange(prop))

## keep cells in the clusters that have more than 10% of OK cells
good_clusters <- t1 %>% filter(prop > 0.10) %>% pull(seurat_clusters)

## export unfiltered per-cell QC metric table
#obj_merged@meta.data = obj_merged[[]] #%>% relocate(dropletQC.keep, .after = 'dropletQC.nucFrac')
#save_qcTale_fn = here('data/tidy_data/tables', 
 #                     paste0("Jarvis_Songbird_unfiltered_QC_table_N",num_samples,'.txt.gz'))

#write_tsv(obj_merged@meta.data, save_qcTale_fn)

## subset cells to those not predicted low QC or doublet
obj_filtered = subset(obj_merged, subset = miQC.keep == "keep" & 
                        scds.keep == "keep" & #& dropletQC.keep == 'cell' &
                        seurat_clusters %in% good_clusters)

## recompute PCA and UMAP embedding post-filtering
obj_filtered = obj_filtered %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(dims = 1:30, verbose = TRUE) %>% 
  FindClusters(resolution = 0.15, algorithm = 2, verbose = TRUE)


#########################################################
# 6) save projects, as h5 object for on-disk computation:
#save_proj_h5_fn = here('data/tidy_data/Seurat_projects', 
#                       paste0("comb_N",num_samples,'.h5Seurat'))
#obj_filtered_test <- Convert(obj_filtered, dest = "h5ad", overwrite = TRUE)
#SaveH5Seurat(obj_filtered, filename = save_proj_h5_fn,overwrite = TRUE)

## save normalized, UMAP embedded, object for downstream analyses
#save_proj_fn = here('data/tidy_data/Seurat_projects', 
#                      paste0("comb_N",num_samples,'.rds'))
saveRDS(obj_filtered, '/projects/pfenninggroup/singleCell/Jarvis_songbird_snRNA-seq/data/tidy_data/final.rds')


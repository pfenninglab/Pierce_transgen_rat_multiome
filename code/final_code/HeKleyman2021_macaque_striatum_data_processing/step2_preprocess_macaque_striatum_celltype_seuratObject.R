## conda activate r4
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

DATADIR='data/tidy_data/HeKleyman2021_macaque_striatum_data_processing'
plan("multicore", workers = 16)
options(future.globals.maxSize = 20 * 1024^3)

###########################################################
## 0) # Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  file = '/home/bnphan/resources/genomes/GRCh38.p13/ENSEMBL_GRCh38.p13_genes_to_Orthologous_mouse_genes.tsv'
  orthologs = read_tsv(file, show_col_types = FALSE) %>% 
    rename_with(make.names) %>% 
    dplyr::filter(Mouse.homology.type== 'ortholog_one2one') %>% 
    dplyr::filter(Gene.name %in% x) %>% 
    dplyr::select(Gene.name, Mouse.gene.name)
  
  # Print the first 6 genes found to the screen
  print(head(orthologs))
  return(orthologs)
}

# RenameGenesSeurat from https://github.com/satijalab/seurat/issues/1049
RenameGenesSeurat <- function(obj, newnames) { 
  # Replace gene names in different slots of a Seurat object. Run this before integration. 
  # It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

#######################################################################
# 1) read in labeled monkey all nuclei dataset, He, Kleyman et al. 2021
## read in the all nuclei object
monkey_all = here(DATADIR, 'rdas', 'GSE167920_Results_full_nuclei_processed_final.rds') %>% readRDS() 
DefaultAssay(monkey_all) <- 'RNA'
monkey_all[['integrated']] <- NULL

## subset to just NAcc which is the region in this study
table(monkey_all$region_name)
monkey_all = monkey_all[,monkey_all$region_name == 'nacc']

## switch out human gene names w/ orthologous mouse gene names, based on ENSEMBL
hgToMm = convertHumanGeneList(rownames(monkey_all)) %>% deframe()
monkey_all = monkey_all[names(hgToMm), ]
monkey_all = RenameGenesSeurat(monkey_all, newnames = hgToMm[rownames(monkey_all)])
head(rownames(monkey_all))
monkey_all

## split by monkey, brain region, and compute SCTransform
monkey_all$split = with(monkey_all[[]], paste(monkey, region_name, gsub('[ACGT-]','', colnames(monkey_all)),sep = '-'))
monkey_all_list <- SplitObject(monkey_all, split.by = "split")
monkey_all_list <- lapply(X = monkey_all_list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = monkey_all_list, nfeatures = 3000)
monkey_all_list <- PrepSCTIntegration(object.list = monkey_all_list, anchor.features = features)
monkey_all_list <- lapply(X = monkey_all_list, FUN = RunPCA, features = features)

## co-embed data into 1 SCTransformed space
monkey_all.anchors <- FindIntegrationAnchors(object.list = monkey_all_list, normalization.method = "SCT", anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
monkey_all.combined.sct <- IntegrateData(anchorset = monkey_all.anchors, normalization.method = "SCT", dims = 1:30)
monkey_all.combined.sct <- RunPCA(monkey_all.combined.sct, verbose = FALSE)
monkey_all.combined.sct <- RunUMAP(monkey_all.combined.sct, reduction = "pca", dims = 1:30)

## save SCTransformed all counts matrix
save_all_h5 = here(DATADIR, 'rdas', 'GSE167920_Results_full_nuclei_processed_final_NAc.mmGenes.h5Seurat')
SaveH5Seurat(monkey_all.combined.sct, filename = save_all_h5, overwrite = TRUE)


#######################################################
# 2) read in labeled monkey MSN subtype dataset, He, Kleyman et al. 2021
## read in the MSN sub-type  object
monkey_msn = here(DATADIR, 'rdas', 'GSE167920_Results_MSNs_processed_final.rds') %>% 
  readRDS()
DefaultAssay(monkey_msn) <- 'RNA'
monkey_msn[['integrated']] <- NULL

## subset to just NAcc which is the region in this study
table(monkey_msn$region_name)
monkey_msn = monkey_msn[,monkey_msn$region_name == 'nacc']

## switch out human gene names w/ orthologous mouse gene names, based on ENSEMBL
hgToMm = convertHumanGeneList(rownames(monkey_msn)) %>% deframe()
monkey_msn = monkey_msn[names(hgToMm), ]
monkey_msn = RenameGenesSeurat(monkey_msn, newnames = hgToMm[rownames(monkey_msn)])
head(rownames(monkey_msn))
monkey_msn

## split by monkey and compute SCTransform
monkey_msn_list <- SplitObject(monkey_msn, split.by = "monkey")
monkey_msn_list <- lapply(X = monkey_msn_list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = monkey_msn_list, nfeatures = 3000)
monkey_msn_list <- PrepSCTIntegration(object.list = monkey_msn_list, anchor.features = features)
monkey_msn_list <- lapply(X = monkey_msn_list, FUN = RunPCA, features = features)

## co-embed data into 1 SCTransformed space
monkey_msn.anchors <- FindIntegrationAnchors(
  object.list = monkey_msn_list, normalization.method = "SCT",
  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
monkey_msn.combined.sct <- IntegrateData(anchorset = monkey_msn.anchors, normalization.method = "SCT", dims = 1:30)
monkey_msn.combined.sct <- RunPCA(monkey_msn.combined.sct, verbose = FALSE)
monkey_msn.combined.sct <- RunUMAP(monkey_msn.combined.sct, reduction = "pca", dims = 1:30)

## save SCTransformed MSN counts matrix
save_msn_h5 = here(DATADIR, 'rdas', 'GSE167920_Results_MSNs_processed_final_NAc.mmGenes.h5Seurat')
SaveH5Seurat(monkey_msn.combined.sct, filename = save_msn_h5, overwrite = TRUE)

save_msn_sce = here(DATADIR, 'rdas', 'GSE167920_Results_MSNs_processed_final_NAc.mmGenes.sce.rds')
monkey_msn.sce = as.SingleCellExperiment(monkey_msn.combined.sct, assay = 'RNA')
saveRDS(monkey_msn.sce, file = save_msn_sce)




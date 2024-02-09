library(scds) ## doublet removal
library(SingleCellExperiment) ## works w/ SCDS
library(SeuratDisk)
library(DropletUtils)
library(Seurat)
library(SoupX)
library(Matrix)
library(here)
library(tidyverse)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)


##################################
# Ambient RNA removal by SoupX
counts_filtered_fn = file.path(list.dirs(countsDIR, recursive = FALSE), 'filtered_feature_bc_matrix')
names(counts_filtered_fn) = ss(counts_raw_fn, '/', 9)

counts_raw_fn = file.path(list.dirs(countsDIR, recursive = FALSE), 'raw_feature_bc_matrix')
names(counts_raw_fn) = ss(counts_raw_fn, '/', 9)

for(sample in names(counts_raw_fn)){
  de_soup_counts_dn = gsub('filtered','SoupX', counts_filtered_fn[sample])
  if(dir.exists(de_soup_counts_dn)){
    print(paste('SoupX already computed:', de_soup_counts_dn))
  } else{
    print(paste('Reading in counts for:', sample,'.'))
    ## read in the filtered counts and raw counts
    ## https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html
    toc = Seurat::Read10X(counts_filtered_fn[sample])[["Gene Expression"]]
    tod = Seurat::Read10X(counts_raw_fn[sample])[["Gene Expression"]]
    
    ## routine Seurat clustering for SoupX, will be redone later
    ## https://satijalab.org/seurat/articles/essential_commands.html
    print(paste('Processing counts for SoupX.'))
    obj <- CreateSeuratObject(counts = toc) %>%
      NormalizeData(verbose = F) %>%
      FindVariableFeatures(verbose = F) %>%
      ScaleData(verbose = F) %>%
      RunPCA(verbose = F) %>%
      FindNeighbors(verbose = F) %>%
      FindClusters(algorithm = 2, resolution = 0.5, verbose = F)
    
    ## give SoupX the filtered counts, raw counts, and highly required basic clustering
    print(paste('Estimating ambient RNA w/ SoupX.'))
    sc = SoupChannel(tod, toc)
    sc = setClusters(sc, setNames(obj$seurat_clusters, colnames(obj)))
    sc = autoEstCont(sc, doPlot = F)
    out = adjustCounts(sc, roundToInt=TRUE)
    
    ## adjust filtered cell counts for ambient RNA
    print(paste('Corrected counts output to:', de_soup_counts_dn))
    
    ## string split (ss function) assumes barcode in the 2 slot of the column name
    write10xCounts(de_soup_counts_dn, out, barcodes = ss(colnames(out),'_',2), 
                   version = '3', overwrite = TRUE)
  }
}

# read in srna-seq matrix files, these folders should have a barcodes.tsv, genes.tsv, and matrix.mtx files
countsDIR= here('data/raw_data/CellRanger_outs')
counts_fn = file.path(list.dirs(countsDIR, recursive = FALSE), 'SoupX_feature_bc_matrix')
names(counts_fn) = ss(counts_fn, '/', 9)
counts_fn

## read in GEM wells and compute scds doublet scores
countList <- counts_fn %>% lapply(Read10X)
# countList <- lapply(countList, function(x) x[["Gene Expression"]])
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
sce_fn = file.path('data/tidy_data/rdas', 'soupx_multiomeRNA_transgen_rat.sce.rds')
saveRDS(sce, sce_fn)

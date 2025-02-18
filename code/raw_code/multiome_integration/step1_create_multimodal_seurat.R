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
library(EnsDb.Mmusculus.v79)
library(BSgenome.Rnorvegicus.UCSC.rn7)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

#########################################################
# 0) Seurat uses the future package for parallelization
plan("multicore", workers = 12)
options(future.rng.onMisuse = 'ignore')
options(future.globals.maxSize = 60 * 1024^3)


###########################################################################
## 1)  import Ensembl database from mouse and lift over to rn7 coordinates
annotations <- import('/home/bnphan/resources/genomes/rn7/rn7_liftoff_mm10_RefSeq.gtf')
seqinfo(annotations) = seqinfo(BSgenome.Rnorvegicus.UCSC.rn7)[seqlevels(annotations)]
annotations <- keepStandardChromosomes(annotations, pruning.mode = "coarse")
mcols(annotations)[,c('score', 'phase', 'source')] = NULL
table(annotations$type)

##############################
## 2) load in the Seurat object
obj = here::here('data/tidy_data/Seurat_projects',
                 "Rat_transgen_multiomeRNA_refined_all_SeuratObj_N5.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')

obj = RenameCells(obj, new.names = ss(colnames(obj), '_', 2))
head(colnames(obj))
colnames(obj) = ss(Cells(obj), '_', 2)

## load in the ArchR project
addArchRThreads(threads = 8)
PROJDIR=here::here("data/tidy_data/ArchRProjects/Rat_Transgen_NAc_multiome")
proj = loadArchRProject(path = PROJDIR)

# 11190 cells have both ATAC and RNA profiles
table( proj$cellNames %in% Cells(obj))
table( Cells(obj) %in% proj$cellNames )

## subset ArchR project to cells in the 
proj <- proj[proj$cellNames %in% colnames(obj) ]
obj = obj[,proj$cellNames]

DefaultAssay(obj) <- "RNA"
obj <- obj %>% 
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

##############################
# 3) create the Chromatin Assay
pkm <- getPeakMatrix(proj)
peaks = getPeakSet(proj)
seqinfo(peaks) = seqinfo(BSgenome.Rnorvegicus.UCSC.rn7)[seqlevels(peaks)]
frag_list = unique(proj$Sample) %>% 
  here('data/raw_data/fragments', . , 'outs/fragments.tsv.gz') %>% 
  lapply(CreateFragmentObject)
         
chrom = CreateChromatinAssay( counts = pkm,
  min.cells = -1, min.features = -1,
  ranges = peaks, fragments = frag_list, genome = seqinfo(peaks),
  annotation = annotations, sep = c("-", "-"),
  validate.fragments = F, verbose = TRUE)

# add the ATAC data to the seurat object and perform dimensionality reduction
head(colnames(chrom))

obj[['ATAC']] = chrom

DefaultAssay(obj) <- "ATAC"
obj <- RunTFIDF(obj) %>% FindTopFeatures(min.cutoff = 'q0') %>% RunSVD() %>% 
   RunUMAP(reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac",
           reduction.key = "atacUMAP_")

# 4) create the motif Assay and perform dimensionality reduction
se_motif <- getMatrixFromProject(proj, useMatrix = "MotifMatrix")
motif_obj <- CreateAssayObject(data = assays(se_motif)$deviations)
head(colnames(motif_obj))

obj[['Motif']] = motif_obj
DefaultAssay(obj) <- "Motif"
VariableFeatures(obj) <- rownames(obj)
obj <- ScaleData(obj) %>% 
  RunPCA(features = rownames(obj), reduction.key = "PCM_", 
         reduction.name = 'pca.motif') %>% 
  RunUMAP(reduction = 'pca.motif', dims = 1:40, reduction.name = "umap.motif",
          reduction.key = "motifUMAP_")

#####################################
## 6) conduct multimodal WNN integration
DefaultAssay(obj) <- "RNA"
obj <- obj %>% 
  FindMultiModalNeighbors(reduction.list = list("pca", "lsi", 'pca.motif'),
                          dims.list = list(1:50, 2:50, 1:40))
obj <- obj %>% 
  RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", 
          reduction.key = "wnnUMAP_")
obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

table(obj$seurat_clusters, obj$cluster_rat)

#####################################
## 7) save the Seurat/Signac object
motimodal_fn =  here::here('data/tidy_data/Seurat_projects',
                           "Rat_transgen_multimodal_SeuratObj_N5.rds")
saveRDS(obj, motimodal_fn)


##########################################################################
## 8) subset to the neurons only and create similar multimodal object
ind_keep = str_detect(obj$cluster_rat, 'MSN|INT')
table(ind_keep, obj$cluster_rat)
obj_neuron = obj[,ind_keep]

# recalculate the variable features for dimensionality reduction
DefaultAssay(obj_neuron) <- "RNA"
obj_neuron <- obj_neuron %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

DefaultAssay(obj_neuron) <- "ATAC"
obj_neuron <- obj_neuron %>% 
  RunTFIDF() %>% FindTopFeatures(min.cutoff = 'q0') %>% RunSVD() %>% 
  RunUMAP(reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac",
          reduction.key = "atacUMAP_")

DefaultAssay(obj_neuron) <- "Motif"
obj_neuron <- obj_neuron %>% 
  RunPCA(features = rownames(obj_neuron), reduction.key = "PCM_", 
         reduction.name = 'pca.motif') %>% 
  RunUMAP(reduction = 'pca.motif', dims = 1:40, reduction.name = "umap.motif",
          reduction.key = "motifUMAP_")

# recompute the WNN integration w/ just the neurons
DefaultAssay(obj_neuron) <- "RNA"
obj_neuron <- obj_neuron %>% 
  FindMultiModalNeighbors(reduction.list = list("pca", "lsi", 'pca.motif'),
                          dims.list = list(1:50, 2:50, 1:40))
obj_neuron <- obj_neuron %>% 
  RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", 
          reduction.key = "wnnUMAP_")
obj_neuron <- FindClusters(obj_neuron, graph.name = "wsnn", algorithm = 3)

#####################################
## 7) save the Seurat/Signac object
motimodal_fn2 =  here::here('data/tidy_data/Seurat_projects',
                           "Rat_transgen_multimodal_SeuratObj_N5.neuron.rds")
saveRDS(obj_neuron, motimodal_fn2)



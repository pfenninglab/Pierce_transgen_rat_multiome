## packages for data table processing 
library(here)
library(tidyverse)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(future)
library(cowplot)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

FIGDIR='figures/exploratory/multiome_integration'
dir.create(here(FIGDIR, 'plots'), recursive = T, showWarnings = F)

######################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 12 cores
plan("multicore", workers = 4)
options(future.globals.maxSize = 100 * 1024^3)
options(future.rng.onMisuse = 'ignore')

celltypes = c('D1-MSN','D1/D3-MSN','D2-MSN','D1-ICj-MSN','D1-NUDAP-MSN','D1/D2-Hybrid-MSN', 
              'Pvalb-INT','Sst-INT', 'Astrocytes', 'Microglia', 'Oligos', 'Oligos_Pre')

cols_rat = ArchR::paletteDiscrete(celltypes)

######################################################
# 1) Seurat uses the future package for parallelization
## load in the full dataset for label refinement
obj = here('data/tidy_data/Seurat_projects', 
           "Rat_transgen_multimodal_SeuratObj_N5.rds") %>% readRDS()
Idents(obj) = 'cluster_rat'
table(Idents(obj) %in%celltypes )
levels(obj) = celltypes

## load in the MSN dataset for label refinement
obj_neuron = here('data/tidy_data/Seurat_projects', 
                       "Rat_transgen_multimodal_SeuratObj_N5.neuron.rds") %>% 
  readRDS()
Idents(obj_neuron) = 'cluster_rat'
levels(obj_neuron) = celltypes


######################################################
# 2) plot the UMAPs for the multimodal data integration
font_size = 8
my_theme = theme(plot.title = element_text(hjust = 0.5, size = font_size + 2),
                   axis.title = element_text(size = font_size), 
                   axis.text = element_text(size = font_size))

p1 <- DimPlot(obj, reduction = "umap.rna", group.by = "cluster_rat", 
              label = T, label.size = 2.5, repel = T, cols = cols_rat) + 
  ggtitle("RNA")+ NoLegend() + my_theme
p2 <- DimPlot(obj, reduction = "umap.atac", group.by = "cluster_rat", 
              label = F, label.size = 2.5, repel = T, cols = cols_rat) + 
  ggtitle("ATAC")+ NoLegend() + my_theme
p3 <- DimPlot(obj, reduction = "umap.motif", group.by = "cluster_rat", 
              label = F, label.size = 2.5, repel = T, cols = cols_rat) + 
  ggtitle("Motif")+ NoLegend() + my_theme
p4 <- DimPlot(obj, reduction = "wnn.umap", group.by = "cluster_rat", 
              label = T, label.size = 2.5, repel = T, cols = cols_rat) + 
  ggtitle("Multi-modal") + NoLegend() + my_theme

pdf(here(FIGDIR, 'plots', 'figure1_multimodal_rat_transgen.all.pdf'),
    width = 8, height = 2.5)
p1 + p2 + p3 + p4 + plot_layout(nrow = 1)
dev.off()



######################################################
# 3) plot the UMAPs for the multimodal data integration neurons only
p5 <- DimPlot(obj_neuron, reduction = "umap.rna", group.by = "cluster_rat", 
              label = T, label.size = 2.5, repel = T, cols = cols_rat) + 
  ggtitle("RNA")+ NoLegend() + my_theme
p6 <- DimPlot(obj_neuron, reduction = "umap.atac", group.by = "cluster_rat", 
              label = F, label.size = 2.5, repel = T, cols = cols_rat) + 
  ggtitle("ATAC")+ NoLegend() + my_theme
p7 <- DimPlot(obj_neuron, reduction = "umap.motif", group.by = "cluster_rat", 
              label = F, label.size = 2.5, repel = T, cols = cols_rat) + 
  ggtitle("Motif")+ NoLegend() + my_theme
p8 <- DimPlot(obj_neuron, reduction = "wnn.umap", group.by = "cluster_rat", 
              label = T, label.size = 2.5, repel = T, cols = cols_rat) + 
  ggtitle("Multi-modal") + NoLegend() + my_theme

pdf(here(FIGDIR, 'plots', 'figure1_multimodal_rat_transgen.neuron.pdf'),
    width = 8, height = 2.5)
p5 + p6 + p7 + p8 + plot_layout(nrow = 1)
dev.off()

#######################################################################
# 4) plot the WNN cell type weights for the multimodal data integration
p9 = VlnPlot(obj, features = "RNA.weight", group.by = 'cluster_rat',
        sort = F, pt.size = 0, cols = cols_rat) + NoLegend() + my_theme + 
  theme(axis.text.x = element_blank()) + 
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'black')
p10 = VlnPlot(obj, features = "ATAC.weight", group.by = 'cluster_rat',
              sort = F, pt.size = 0, cols = cols_rat) + NoLegend() + my_theme+ 
  theme(axis.text.x = element_blank())+ 
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'black')
p11 = VlnPlot(obj, features = "Motif.weight", group.by = 'cluster_rat',
              sort = F, pt.size = 0, cols = cols_rat) + NoLegend() + my_theme+ 
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'black')

pdf(here(FIGDIR, 'plots', 'sfig1_multimodal_rat_transgen.WNNweights.pdf'),
    width = 8, height = 6)
p9 + p10 + p11 + plot_layout(ncol = 1)
dev.off()



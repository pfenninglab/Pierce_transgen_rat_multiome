## packages for data table processing 
library(ArchR)
library(tidyverse)
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

#########################################
## 1) import the metadata with the labels
save_meta_fn = here('data/tidy_data/tables',"Rat_transgen_multiomeRNA_refined_all_SeuratObj_N5.txt.gz")
df_meta = read_tsv(save_meta_fn, show_col_types = FALSE) %>% 
  dplyr::select(-c(nCount_RNA:Sample, Sire:celltype2)) %>% 
  column_to_rownames('Barcode')

with(df_meta, table(cluster_rat, cluster_macaque))

#####################################################
## 2) load the ATAC data to add the cell cluster types
proj = here::here('data/tidy_data/ArchRProjects','Rat_Transgen_NAc_scATAC') %>% 
  loadArchRProject(ARCHDIR)

proj2 = here::here('data/tidy_data/ArchRProjects','Rat_Transgen_NAc_multiome') %>% 
  loadArchRProject()

df_meta = df_meta[rownames(df_meta)%in% getCellNames(proj), ]
df_meta2 = df_meta[rownames(df_meta)%in% getCellNames(proj2), ]

## add the cell type labels to the ArchR project
for(col in names(df_meta)){
  ## for the scATAC-seq
  proj = addCellColData(proj, data = df_meta[,col], force = TRUE,
                        name = col, cells = rownames(df_meta))
  
  ## for the multiome
  proj2 = addCellColData(proj2, data = df_meta2[,col], force = TRUE, 
                        name = col, cells = rownames(df_meta2))
}

## check the clusters
table(proj$Sample, proj$cluster_macaque)
table(proj$Sample, proj$cluster_rat)

## save project
proj = saveArchRProject(proj)



#####################################################
## 2) load the ATAC data to add the cell cluster types

df_meta = df_meta[rownames(df_meta)%in% getCellNames(proj), ]

## add the cell type labels to the ArchR project
for(col in names(df_meta)){
  proj = addCellColData(proj, data = df_meta[,col],
                        name = col, cells = rownames(df_meta), force = TRUE)
}

## check the clusters
table(proj$Sample, proj$cluster_macaque)
table(proj$Sample, proj$cluster_rat)

## save project
proj = saveArchRProject(proj)

## packages for data table processing 
library(here)
library(tidyverse)
library(Seurat)
library(Signac)
library(future)
library(presto)
library(metap)
library(multtest)
library(qvalue)

DATADIR='data/tidy_data/tables'
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

plan("multicore", workers = 20)
options(future.globals.maxSize = 25 * 1024^3)

######################################################
# 1) Seurat uses the future package for parallelization
## load in the full dataset for label refinement
obj = here('data/tidy_data/Seurat_projects', 
           "Rat_transgen_multimodal_SeuratObj_N5.rds") %>% readRDS()
DefaultAssay(obj) = 'peak'
Idents(obj) = 'Sire'

# Find differential ATAC peaks between conditions
df_cocVsal = FindConservedMarkers( 
  obj, ident.1 = 'Cocaine', ident.2 = 'Saline', 
  grouping.var = 'cluster_rat', assay = "ATAC", slot = "data", 
  min.cells.group = 10, min.cells.feature = 5, 
  test.use = "wilcox", densify = F, min.pct=0, logfc.threshold = 0, verbose = T)

df_methVsal = FindConservedMarkers( 
  obj, ident.1 = 'Methamphetamine', ident.2 = 'Saline', 
  grouping.var = 'cluster_rat', assay = "ATAC", slot = "data", 
  min.cells.group = 10, min.cells.feature = 5,
  test.use = "wilcox", densify = T, min.pct=0, logfc.threshold = 0, verbose = T)

df_methVcoc = FindConservedMarkers( 
  obj, ident.1 = 'Methamphetamine', ident.2 = 'Cocaine', 
  grouping.var = 'cluster_rat', assay = "ATAC", slot = "data", 
  min.cells.group = 10, min.cells.feature = 5,
  test.use = "wilcox", densify = T, min.pct=0, logfc.threshold = 0, verbose = T)


save(df_cocVsal, df_methVsal, df_methVcoc, 
     file = here('data/tidy_data/differential_expression/rdas/diff_peak_3-way.rda'))

#######################################################################################
# 2) Combine the dataframes of differential peaks

## Helper variables for reshaping dataframe
celltypeMatch = obj$cluster_rat %>% unique() %>% paste0('^',.,'_') %>% paste(collapse = '|')
metricMatch = names(df_methVcoc) %>% grep(celltypeMatch, ., value = T) %>%  
  str_replace(celltypeMatch, '') %>% unique() %>% paste0('_', ., '$') %>% paste(collapse = '|')

## ATAC average activity of each ATAC in each cell type
avg_activity = AverageExpression(obj, assay = "ATAC", group.by = "cluster_rat")[[1]]
avg_activity = avg_activity %>% as.data.frame() %>% rownames_to_column('peak') %>% 
  pivot_longer(-peak, names_to = 'celltype', values_to = 'AveActivity')

## Reshape the differential ATAC table
df = list('Coc_vs_Sal' = df_cocVsal,
          'Met_vs_Sal' = df_methVsal,
          'Met_vs_Coc' = df_methVcoc) %>% 
  lapply(rownames_to_column, 'peak') %>% bind_rows(.id = 'group') %>% 
  pivot_longer(cols = -c(group, peak)) %>% 
  mutate(metric = str_replace_all(name, celltypeMatch, ''),
         celltype = str_replace_all(name, metricMatch, '')) %>% 
  dplyr::select(-name) %>% 
  filter(!metric %in% c('max_pval', 'minimump_p_val')) %>% 
  pivot_wider(values_from = 'value', names_from = 'metric')

## Add average activity and adjust p-values
df = df %>% inner_join(avg_activity) %>% 
  group_by(group, celltype) %>% 
  mutate(p_val_adj = qvalue(p_val)$qvalues,
         p_val_bonf = p.adjust(p_val, 'bonferroni')) %>% 
  ungroup()

table(df$p_val_adj < 0.05, df$group) 
table(df$p_val_adj < 0.05, df$group, df$celltype) 

################################################################
# 3) Analyze differential peaks 
df_wide = df %>% dplyr::select(-c(pct.1, pct.2)) %>% 
  relocate(c(ATAC, celltype, AveActivity), .before = everything()) %>% 
  pivot_wider(id_cols = ATAC:AveActivity, names_from = 'group', 
              values_from = c(p_val:p_val_bonf)) %>% 
  mutate(
    ## Calculate effect sizes
    effect_Met = -log10(p_val_Met_vs_Sal) * sign(avg_log2FC_Met_vs_Sal),
    effect_Coc = -log10(p_val_Coc_vs_Sal) * sign(avg_log2FC_Coc_vs_Sal),
    ## Calculate interaction score
    interaction_MvC = (effect_Met - effect_Coc)) %>% 
  arrange(desc(abs(interaction_MvC)))

## List of dataframes by cell type
df_list = df_wide %>% split(f = .$celltype)
names(df_list) = names(df_list) %>% make.names()

## Find significant differential peaks
alpha = 0.05
df_interact = df_wide %>% 
  filter(p_val_adj_Coc_vs_Sal < alpha | 
           p_val_adj_Met_vs_Sal < alpha & 
           p_val_adj_Met_vs_Coc < alpha) %>% 
  ## the interaction meth vs. coc must be significant
filter(abs(interaction_MvC)>1 )


######################################################
# 4) Summarize results and export
## Count differential peaks
df %>% group_by(group, celltype) %>% 
  summarise(numDM = sum(p_val_adj < alpha & abs(avg_log2FC) > 0)) %>% 
  pivot_wider(names_from = 'group', values_from = 'numDM')

df %>% group_by(group, celltype) %>% 
  summarise(numDM = sum(p_val_bonf < alpha & abs(avg_log2FC) > 0)) %>% 
  pivot_wider(names_from = 'group', values_from = 'numDM')

## Summarize interactions
table(df_interact$celltype)
df_interact %>% filter(abs(interaction_MvC) > 1) %>% 
  group_by(celltype) %>% summarise(numDM_interact = n())

## Save results
dir.create(here('data/tidy_data/differential_expression/rdas'), showWarnings = F)
saveRDS(df_wide, file = here('data/tidy_data/differential_expression/rdas/diff_peak_3-way.rds'))

dir.create(here('data/tidy_data/differential_expression/tables'), showWarnings = F)
df_wide %>% filter(celltype == 'D1-ICj-MSN') %>% 
  writexl::write_xlsx(here('data/tidy_data/differential_expression/tables',
                                     'differential_peaks_3-way_sires.D1-ICj-MSN.xlsx'))
df_interact %>% writexl::write_xlsx(here('data/tidy_data/differential_expression/tables',
                                         'differential_peaks_top_interaction_MvC.xlsx'))




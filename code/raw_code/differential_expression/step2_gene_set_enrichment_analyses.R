## packages for data table processing 
library(here)
library(tidyverse)
library(data.table)

library(fgsea)
library(swfdr)
library(msigdbr)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

#######################################################
# 0) Seurat uses the future package for parallelization
plan("multicore", workers = 28)
options(future.globals.maxSize = 100 * 1024^3)
options(future.rng.onMisuse = 'ignore')

#################################################################################
# 1) read in the DEG lists and rank order the genes by each differential comparisons
df_wide = readRDS(here('data/tidy_data/differential_expdf_longsion/rdas/diff_gene_3-way_countsplit.rds'))

metric_match = "p_val_|avg_log2FC_|p_val_adj_|p_val_bonf_"
group_match = "_Met_vs_Sal|_Met_vs_Coc|_Coc_vs_Sal"

## grab the p-values and the log2FC differences
df_long = df_wide %>% pivot_longer(-c(gene:AveExpr)) %>% 
  filter(!grepl('effect|inter',name)) %>% 
  mutate(group = gsub(metric_match, '', name), 
         metric = gsub(group_match, '', name)) %>% 
  dplyr::select(-name) %>% 
  pivot_wider(names_from = 'metric') 

## differentially expressed genes
alpha = 0.001
df_long = df_long %>% filter(p_val_adj < alpha) %>% arrange(p_val)
table( df_long$celltype, df_long$group)

# 2) read in the clustered pathway enrichment tables
curenrich_clustered = 
  here(FIGDIR, 'tables','clustered_gsea_pathway_network_interaction_sets.xlsx') %>% 
  readxl::read_xlsx() %>% mutate(
    leadingEdge = str_split(leadingEdge, ',')
  ) %>% unnest(leadingEdge) %>% dplyr::rename('gene'= 'leadingEdge') %>% 
  dplyr::rename('group' = 'comparison')

########################################
## D1H and Microglia Synapse cluster
curenrich_clustered %>%
  filter(cluster_number %in% c(1)) %>%
  count(gene, celltype, group) %>% 
  arrange(gene) %>%
  inner_join(df_long) %>% 
  filter(n > 1) %>% 
  as.data.frame()

## add description
curenrich_clustered %>%
  filter(cluster_number %in% c(2)) %>%
  count(gene, celltype, group) %>% 
  arrange(gene) %>%
  inner_join(df_long) %>% 
  filter(n > 1, celltype =='Microglia') %>% 
  as.data.frame()


curenrich_clustered %>%
  filter(cluster_number %in% c(2)) %>%
  count(gene, celltype, group) %>% 
  arrange(gene) %>%
  inner_join(df_long) %>% 
  filter(n > 1, celltype =='Oligos_Pre') %>% 
  as.data.frame()


## microglia cluster
curenrich_clustered %>% filter(celltype == 'Microglia')  %>% 
  distinct(celltype, gene) %>% 
  inner_join(df_long) %>% arrange(adj.P.Val.Between)

df_long %>% filter(celltype == 'Microglia') %>% filter(gene %in% c('ADGRB3'))
df_long %>% filter(gene %in% c('NFKBIA'))
df_long %>% filter(gene %in% c('APOE'))

## endothelial cluster
curenrich_clustered %>% filter(cluster_number %in% c(1, 2, 4, 5, 6)) %>%
  count(gene) %>% arrange(desc(n)) %>% filter(n>5) %>% 
  mutate(celltype = 'Endothelial') %>% inner_join(df_long)


## the down-regulated cluster
curenrich_clustered %>% 
  count(gene) %>% arrange(desc(n)) %>% filter(n>4) %>% inner_join(df_long) %>% 
  filter(!celltype %in% c('All', 'Glia')) %>% filter(grepl('GAB|GRM|GRIA|GRIN', gene))

curenrich_clustered %>% filter(cluster_number == 9)  %>% 
  count(gene) %>% arrange(desc(n)) %>% filter(n>=3) %>% pull(gene)

curenrich_clustered %>% filter(cluster_number == 9)  %>% 
  count(gene) %>% arrange(desc(n)) %>% filter(n>=3) %>% inner_join(df_long) %>% 
  filter(!celltype %in% c('All', 'Glia'))


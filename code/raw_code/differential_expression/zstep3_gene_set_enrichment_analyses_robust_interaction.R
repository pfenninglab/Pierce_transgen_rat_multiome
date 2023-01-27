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
df_wide = readRDS(here('data/tidy_data/differential_expression/rdas/diff_gene_3-way_countsplit.rds'))

metric_match = "p_val_|avg_log2FC_|p_val_adj_|p_val_bonf_"
group_match = "_Met_vs_Sal|_Met_vs_Coc|_Coc_vs_Sal"

## grab the p-values and the log2FC differences
deg_rank_list1 = df_wide %>% pivot_longer(-c(gene:AveExpr)) %>% 
  filter(!grepl('effect|inter',name)) %>% 
  mutate(group = gsub(metric_match, '', name), 
         metric = gsub(group_match, '', name)) %>% 
  dplyr::select(-name) %>% 
  pivot_wider(names_from = 'metric') %>% 
  mutate(value = -log10(p_val) * sign(avg_log2FC), 
         group2 = paste(celltype, group, sep = '#')) %>% 
  split(f = .$group2) %>% 
  lapply(function(x){
    ## some of the unexpressed genes
    x %>% filter(value!=0) %>% 
      ## the ranks can't be ties
      mutate(value = jitter(value, factor = 1e-2)) %>% 
      dplyr::select(gene, value) %>% arrange(value) %>% deframe() 
  })  

## grab the ranks by the interaction terms
deg_rank_list2 = df_wide %>% 
  ## use the interaction metric between meth vs. saline VS. coc vs. saline
  mutate(value = interaction_MvC, group2 = paste0(celltype, '#interaction')) %>% 
  split(f = .$group2) %>% 
  lapply(function(x){
    ## some of the unexpressed genes
   x %>% filter(value!=0) %>% 
      ## the ranks can't be ties
      mutate(value = jitter(value, factor = 1e-2)) %>% 
   dplyr::select(gene, value) %>% arrange(value) %>% deframe() 
  })

deg_rank_list = c(deg_rank_list1, deg_rank_list2)
lengths(deg_rank_list) # 36


#############################################
## 2) get gene ontologies, use human genes
## grab the H, Hallmark set of gene pathways
## grab the C2, which is the curated canonical pathway sets
## grab the C5, which is the Gene Ontology sets
## https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
pathways_df =  bind_rows(msigdbr("mouse", category="H"), 
                         msigdbr("mouse", category="C2"), 
                         msigdbr("mouse", category="C5"))

## get the SynGO gene ontologies, use human genes
syngo_df = readRDS(here('data/tidy_data/SynGO_bulk_download_release_20210225',
                        'rdas', 'syngo_annotations.rds'))
pathways_df = rbindlist(list(pathways_df, syngo_df), fill = T) 

## reshape/label for pathway naming purposes
pathways <-pathways_df %>% 
  mutate(
    gs_subcat = ifelse(is.na(gs_subcat) | gs_subcat == '', gs_cat, gs_subcat),
    gs_name = paste(gs_subcat, gs_name, sep ='#')) %>% 
  split(x = .$gene_symbol, f = .$gs_name)

## exclude the really really big gene sets
lengths(pathways) %>% summary()
pathways = pathways[lengths(pathways)<500]
length(pathways) # 21437

table(pathways_df$gs_cat)

pathways_df2 = pathways_df %>% 
  dplyr::select(gs_subcat, gs_name, gs_description) %>% distinct() %>% 
  dplyr::rename('pathway_group' = 'gs_subcat', 'pathway' = 'gs_name',
                'description' = 'gs_description')

#################################################################
## 3) conduct the GSEA analyses across the different comparisons
gsea_list = lapply(deg_rank_list, fgsea, pathways = pathways,
                   minSize=15, ## minimum gene set size
                   maxSize=400) ## maximum gene set size

alpha = 0.05
gsea_df = gsea_list %>% rbindlist(idcol = 'group') %>% 
  arrange(pval) %>% filter(!is.na(pval)) %>% 
  mutate(
    MSigDb_Group = ss(pathway, '#', 1), 
    pathway = ss(pathway, '#', 2),
    padj = lm_qvalue(pval, X=size)$q, 
    celltype = group %>% ss('#', 1), 
    comparison = group %>% ss('#', 2), 
    leadingEdge = map_chr(leadingEdge, paste, collapse = ',')) %>% 
  ## use a lower cutoff of for replication
  inner_join(pathways_df2) %>% filter(pval < alpha) %>% 
  dplyr::select(-group) %>% 
  relocate(comparison, celltype, MSigDb_Group, pathway, description, .before= everything()) 

## 
gsea_df %>% split(f = .$comparison) %>% 
  writexl::write_xlsx( here('data/tidy_data/differential_expression', 
                            'tables/GSEA_enrichment_msigdb_H_C2_C5_SynGO.xlsx') )
 

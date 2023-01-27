## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RColorBrewer)
library(rcartocolor)
library(ggpubr)

library(data.table)
library(fgsea)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)


## make for this subdirs
PLOTDIR='figures/exploratory/Piechota_et_al_mouse_striatum_on_drugs'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)


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


## read in the gene sets from the Piechota et al. paper
piechota_list = here('data/tidy_data/Piechota_et_al_mouse_striatum_on_drugs/tables',
               "13059_2010_2341_MOESM2_ESM.XLS") %>%
  readxl::read_xls('patterns_gene_expression_values') %>% 
  rename_with(make.names) %>% dplyr::rename('Mouse.gene.name' = 'Gene.Symbol') %>% 
  dplyr::select(Mouse.gene.name:B3_extended) %>% 
  pivot_longer(cols = A:B3_extended,names_to = 'pattern', values_to = 'value' ) %>% 
  filter(!is.na(value) & grepl('ext', pattern)) %>% 
  mutate(pattern = ss(pattern, '_') %>% paste0('Piechota2010.',.)) %>% 
  split(x = .$Mouse.gene.name, f = .$pattern)
    

## conduct the GSEA analyses
gsea_list = lapply(deg_rank_list, fgsea, pathways = piechota_list,
                   minSize=5, ## minimum gene set size
                   maxSize=400) ## maximum gene set size

gsea_df = gsea_list %>% rbindlist(idcol = 'group') %>% arrange(pval) %>% 
  mutate(padj = p.adjust(pval, 'fdr'), 
         celltype = group %>% ss('#', 1), 
         OUD.v.CTL.in = group %>% ss('#', 2), 
         leadingEdge = map_chr(leadingEdge, paste, collapse = ',')) %>% 
  dplyr::select(-group) %>% relocate(celltype, OUD.v.CTL.in, .before= everything())

out_fn = here(PLOTDIR, 'tables', 'GSEA_enrichment_Piechota_et_al_mouse_striatum_on_drugs.xlsx')
gsea_df %>% writexl::write_xlsx(out_fn)


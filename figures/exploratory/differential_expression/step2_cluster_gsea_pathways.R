## packages for data table processing 
## this script was adapted from 
## https://github.com/kowaae22/ClarkLabDocumentation/blob/main/MakeigraphFigures.R
## packages for data table processing 
library(tidyverse)
library(igraph)
library(network)
library(sna)
library(stringr)
library(RColorBrewer)
library(data.table)
library(here)
library(leiden)

"%ni%" <- function(x, table) !(match(x, table, nomatch = 0) > 0)
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
source(here('code/final_code/Rutils/igraph_pathway_clustering.R'))

FIGDIR='figures/exploratory/differential_expression'
dir.create(here(FIGDIR, 'plots', 'gsea'), recursive = T, showWarnings = F)
dir.create(here(FIGDIR, 'tables'), recursive = T, showWarnings = F)

###########################################
# 1) read in the pathway enrichment tables
xl_path = here('data/tidy_data/differential_expression', 
               'tables/GSEA_enrichment_msigdb_H_C2_C5_SynGO.xlsx')
readxl::excel_sheets(path = xl_path)

## read in the comparison of OUD vs. CTL averaged across all subtypes
interaction_df = bind_rows(readxl::read_excel(path = xl_path, sheet = 'Met_vs_Coc'), 
                           readxl::read_excel(path = xl_path, sheet = 'interaction')) %>% filter(padj < 0.05)
gsea_coca = readxl::read_excel(path = xl_path, sheet = 'Coc_vs_Sal')%>% filter(pbon < 0.05)
gsea_meth = readxl::read_excel(path = xl_path, sheet = 'Met_vs_Sal')%>% filter(pbon < 0.05)

## filter out the pathways are in the same direction
gsea_coca2 = gsea_coca %>% 
  mutate(NES_other = gsea_meth$NES[match(pathway, gsea_meth$pathway)]) %>% 
  filter(pathway %ni% gsea_meth$pathway | sign(NES_other) != sign(NES)) %>% 
  dplyr::select(-NES_other) %>% filter(pbon < 0.05) 

gsea_meth2 = gsea_meth %>% 
  mutate(NES_other = gsea_coca$NES[match(pathway, gsea_coca$pathway)]) %>% 
  filter(pathway %ni% gsea_coca$pathway | sign(NES_other) != sign(NES)) %>% 
  dplyr::select(-NES_other) %>% filter(pbon < 0.05) 

## subset by female and male NES to create two pathway networks
with(gsea_coca2, table(celltype, comparison))
net_coca=make_igraph_from_pathways(gsea_coca2)
V(net_coca)$shape = ifelse(V(net_coca)$pathway %in% interaction_df$pathway, 'square', 'circle')

with(gsea_meth2, table(celltype, comparison))
net_meth=make_igraph_from_pathways(gsea_meth2)
V(net_meth)$shape = ifelse(V(net_meth)$pathway %in% interaction_df$pathway, 'square', 'circle')


##################################################
# 2)  add the colors to the NES of either pathways
nes_max = max(c(V(net_coca)$NES, V(net_meth)$NES)) 
nes_min = min(c(V(net_coca)$NES, V(net_meth)$NES))
nes_range = nes_max - nes_min

col = circlize::colorRamp2(seq(nes_min, nes_max, nes_range/10), brewer.pal(11,"PiYG"))
V(net_coca)$color= col(V(net_coca)$NES)
V(net_meth)$color= col(V(net_meth)$NES)


##############################################################################
## trim and cluster the female network; CHANGE target_degree if you have too 
## many edges, clustering will not work + you will get a hairball network
net_coca.sp = trim_edges(net_coca, target_degree = .1)
clp_coca <- cluster_leiden(net_coca.sp, resolution_parameter = 1, 
                             objective_function = 'modularity', n_iterations = 10)
pdf(here(FIGDIR, 'plots','gsea','clustered_coca_gsea_pathway_network.pdf'))
plot(clp_coca, net_coca.sp, vertex.label=NA)
dev.off()

net_coca_list = trim_clustered_nodes(net_coca.sp, clp_coca, min_nodes = 3)
V(net_coca_list$net)$hub_score = hub_score(net_coca_list$net)$vector
V(net_coca_list$net)$label.cex = .5

##############################################################################
## trim and cluster the male network; CHANGE target_degree if you have too 
## many edges, clustering will not work + you will get a hairball network
net_meth.sp = trim_edges(net_meth, target_degree = .1)
clp_meth <- cluster_leiden(net_meth.sp, resolution_parameter = 1, 
                           objective_function = 'modularity',  n_iterations = 10)
pdf(here(FIGDIR, 'plots','gsea','clustered_meth_gsea_pathway_network.pdf'))
plot(clp_meth, net_meth.sp, vertex.label=NA)
dev.off()

net_meth_list = trim_clustered_nodes(net_meth.sp, clp_meth, min_nodes = 3)
V(net_meth_list$net)$mem = V(net_meth_list$net)$mem + max(V(net_coca_list$net)$mem)
V(net_meth_list$net)$hub_score = hub_score(net_meth_list$net)$vector
V(net_meth_list$net)$label.cex = .5


###################################################################
## clean up the labels, merging the pathways from both clusterings
vert.meta =  bind_rows(vertex_attr(net_coca_list$net), 
                       vertex_attr(net_meth_list$net)) %>% as.data.frame()
to_label = vert.meta %>% clean_pathways()
to_label_num = vert.meta %>% dplyr::select(name, mem) %>% deframe()

## make sets of 3 plots for the current female network
in2mm<-25.4
prefix_coca = here(FIGDIR, 'plots','gsea','clustered_coca_gsea_pathway_network')
plot_clusterings_3set(net_coca_list, prefix_coca, to_label, to_label_num, 
                      height = 50/in2mm, width = 50/in2mm)
plot_pt_legend(V(net_meth_list$net)$size, prefix_coca, height = 20/in2mm, width = 10/in2mm)
plot_col_legend(prefix = prefix_coca, col = col, fontsize = 5,
                height = 8/in2mm, width = 30/in2mm)


## make sets of 3 plots for the current male network
prefix_meth = here(FIGDIR, 'plots','gsea','clustered_meth_gsea_pathway_network')
plot_clusterings_3set(net_meth_list, prefix_meth, to_label, to_label_num, 
                      height = 50/in2mm, width = 50/in2mm)

## add the clustered pathways back into the 
curenrich_clustered = bind_rows(gsea_coca2, gsea_meth2, 
                                interaction_df %>% filter(pbon < 0.05) ) %>% 
  mutate(cluster_number = to_label_num[ paste0(celltype, '#', pathway)]) %>% 
  arrange(cluster_number, padj) %>%
  relocate(cluster_number, pathway, description,.after = 'celltype') %>% 
  writexl::write_xlsx(here(FIGDIR, 'tables','clustered_gsea_pathway_network_interaction_sets.xlsx'))

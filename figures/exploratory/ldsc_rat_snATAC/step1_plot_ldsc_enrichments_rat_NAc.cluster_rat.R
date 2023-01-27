#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggsci))
suppressMessages(library(ArchR))
library(wesanderson)
library(ggrepel)
library(ggpubr)
library(here)
# to be run in the root github directory
LABEL='ldsc_rat_snATAC'
DATADIR='data/tidy_data/ldsc_rat_snATAC'

##########################
# read in the GWAS traits
load(here(DATADIR,'rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'), 
         trait = factor(trait, rev(unique(trait))))

celltypes = c('D1.MSN','D1.D3.MSN','D2.MSN','D1.ICj.MSN','D1.NUDAP.MSN','D1.D2.Hybrid.MSN', 
              'Pvalb.INT','Sst.INT', 'Astro', 'Microglia', 'Oligo', 'Oligos_Pre')

col_celltype = c(RColorBrewer::brewer.pal(6, 'Paired'),
                 RColorBrewer::brewer.pal(6, 'Dark2'))
names(col_celltype) = celltypes

#########################################
# read in the LDSC partitioned heritability estimation
enrich_fn =here(DATADIR,'enrichments') %>% 
  list.files(path = ., pattern = '.cell_type_results.txt', full.names = T)
names(enrich_fn) = ss(basename(enrich_fn), '.cell_type_results.txt')
input = lapply(enrich_fn, read_tsv, show_col_types = FALSE) %>% 
  bind_rows(.id = 'file') %>% 
  mutate(celltype = str_replace(Name, 'Rat_Transgen_NAc.', '')) %>% 
  filter(celltype %in% celltypes)
input %>% data.frame() %>% head()

#########################################
## format groupings and calculate conditional cell type enrichment p-value
enrichments = input %>% 
  mutate(
    file = gsub('Rat_Transgen_NAc.', '', file), 
    match = ss(file, '\\.', 1), 
    celltype = gsub('Rat_Transgen_NAc.', '', Name),
    celltype = factor(celltype,celltypes), 
    cell_group = case_when(
      celltype %in% celltypes[1:3] ~ 'Canonical MSN', 
      celltype %in% celltypes[4:6] ~ 'Non-canonical MSN', 
      TRUE ~ 'Interneuron and Glia'
    ), 
    cell_group = factor(cell_group, c('Canonical MSN','Non-canonical MSN', 
                                      'Interneuron and Glia'))) %>%
  inner_join(x = pheno, by = 'match') %>% filter(!grepl('UNK',celltype))

enrichments %>% data.frame() %>% head()

## normalize the coefficients by per SNP heritability
enrichments = enrichments %>%
  mutate(Pbonf = p.adjust(Coefficient_P_value, 'bonferroni'), 
         logPbonf = -log10(Pbonf), 
         FDR = p.adjust(Coefficient_P_value, 'fdr'), 
         logFDR = -log10(FDR), 
         p.signif = Pbonf < 0.05,
         p.signifFDR = FDR < 0.05,
         Coef_norm = Coefficient / h2_perSNP, 
         Coef_norm_se = Coefficient_std_error / h2_perSNP) %>% 
  ungroup()

# look at top signif celltypes
enrichments %>% filter(p.signif) %>% select(match, celltype, Pbonf) %>%
  group_by(match) %>% top_n(1, -Pbonf) %>% ungroup() %>% 
  data.frame() %>% arrange(Pbonf) 


#################################
## make plots for presentation ##
alpha = 0.05
height_ppt = 4.5; width_ppt = 8
height_fig = 5; width_fig = 2.25; font_fig = 5
pal <- wes_palette("Rushmore1", 10, type = "continuous")

# make plots
FIGDIR=here('figures/exploratory/ldsc_rat_snATAC')
dir.create(FIGDIR, showWarnings = F, recursive = T)
plot_fn = here(FIGDIR, 'plots','Rat_Transgen_NAc_ldsc.cluster_rat.bonf.ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn)
for(plotme in list(c('Psych',"SUD", 'Neuro','Sleep'), 
                   c('SU','Degen','Metabolism','Other'))){
  pp = ggplot(data = enrichments %>% filter(group %in% plotme),
              aes(y = trait, x = celltype, fill = logPbonf,color = p.signif)) +
    geom_tile(size = .4) +
    scale_color_manual(values = c('black','#00000015'), breaks = c(TRUE,FALSE),
                       name = paste('P_bonf <',alpha)) + 
    scale_fill_gradientn(colours = pal, name = "-Log10(P bonf)", 
                         limits = c(0, max(enrichments$logPbonf))) +
    facet_grid(group ~ cell_group, scales = 'free', space = 'free') +  
    xlab('Cell type') + ylab('GWAS Trait') + 
    theme_bw(base_size = 9) + 
    guides(colour = guide_legend(override.aes = list(fill = 'white'))) + 
    theme(legend.position = "right", 
          legend.text=element_text(size=8),
          legend.title=element_text(size=9))+ 
    theme(axis.text.x = element_text(angle = -30, hjust = 0))
  print(pp)
}
dev.off()



plot_fn = here(FIGDIR, 'plots','Rat_Transgen_NAc_ldsc.cluster_rat.FDR.ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn)
for(plotme in list(c('Psych',"SUD", 'Neuro','Sleep'), 
                   c('SU','Degen','Metabolism','Other'))){
  pp = ggplot(data = enrichments %>% filter(group %in% plotme),
              aes(y = trait, x = celltype, fill = logFDR,color = p.signifFDR)) +
    geom_tile(size = .4) +
    scale_color_manual(values = c('black','#00000015'), breaks = c(TRUE,FALSE),
                       name = paste('FDR <',alpha)) + 
    scale_fill_gradientn(colours = pal, name = "-Log10(FDR)", 
                         limits = c(0, max(enrichments$logPbonf))) +
    facet_grid(group ~ cell_group, scales = 'free', space = 'free') +  
    xlab('Cell type') + ylab('GWAS Trait') + 
    theme_bw(base_size = 9) + 
    guides(colour = guide_legend(override.aes = list(fill = 'white'))) + 
    theme(legend.position = "right", 
          legend.text=element_text(size=8),
          legend.title=element_text(size=9))+ 
    theme(axis.text.x = element_text(angle = -30, hjust = 0))
  print(pp)
}
dev.off()


#################################
## explort gwas enrichment table ##
dir.create(here('figures/exploratory/ldsc_rat_snATAC', 'tables'), recursive = F, showWarnings = F)
out_fn = here('figures/exploratory/ldsc_rat_snATAC', 'tables', 
              'Rat_Transgen_NAc.GWAS_enrichment_results.cluster_rat.xlsx')
enrichments %>% writexl::write_xlsx(out_fn)



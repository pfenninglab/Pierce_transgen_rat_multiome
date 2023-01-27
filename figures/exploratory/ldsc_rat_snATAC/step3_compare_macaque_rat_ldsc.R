library(tidyverse)
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

# to be run in the root github directory
LABEL='ldsc_rat_snATAC'
DATADIR='data/tidy_data/ldsc_rat_snATAC'

celltypes = c('D1.Shell', 'D2.Shell','D1.ICj','D1.NUDAP','D1.D2.Hybrid', 
              'Interneurons', 'Astro', 'Microglia', 'Oligo', 'Oligos_Pre')

col_celltype = c(RColorBrewer::brewer.pal(10, 'Paired'))
names(col_celltype) = celltypes

########################################
## 1) read in the GWAS enrichments
macaque_df = '/projects/pfenninggroup/singleCell/Macaque_snATAC-seq/macaque_snATAC-seq/Macaque_Multiome_Striatum/figures/exploratory/ldsc_multiomeATAC_Striatum/tables/Stauffer_striatum.GWAS_enrichment_results.xlsx' %>% 
  readxl::read_xlsx() %>% mutate(Species = 'Macaca_mulatta') %>% 
  mutate(celltype = ifelse(celltype =='OPC', 'Oligos_Pre', celltype), 
         ifelse(celltype =='OPC', 'Oligos_Pre', celltype)) %>% 
  filter(celltype %in% celltypes)  %>% 
  dplyr::select(c(trait, celltype, Coefficient:Species))
table(macaque_df$celltype)


rat_df = here('figures/exploratory/ldsc_rat_snATAC', 'tables', 
              'Rat_Transgen_NAc.GWAS_enrichment_results.cluster_macaque.xlsx') %>% 
  readxl::read_xlsx() %>% mutate(Species = 'Rattus_norvegicus') %>% 
  relocate(celltype:cell_group, .after = 'label') %>% 
  dplyr::select(c(trait, celltype, Coefficient:Species))
table(rat_df$celltype)

alpha = 0.05

df_wide = bind_rows(rat_df, macaque_df) %>%
  pivot_wider(id_cols= c(trait, celltype), names_from = Species, 
              values_from = c(Coefficient:Coef_norm_se), values_fn = mean) %>% 
  filter(FDR_Rattus_norvegicus < alpha | FDR_Macaca_mulatta < alpha) %>%
  filter(!is.na(Coef_norm_Macaca_mulatta) & !is.na(Coef_norm_Rattus_norvegicus))
table(df_wide$celltype)

# make plots
height_ppt = 4.5; width_ppt = 8
height_fig = 5; width_fig = 2.25; font_fig = 5

FIGDIR=here('figures/exploratory/ldsc_rat_snATAC')
plot_fn = here(FIGDIR, 'plots','Rat_Transgen_NAc.vs.Stauffer_Monkey_striatum.LDSC.ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn)

pp = ggplot(data = df_wide,
            aes(x = Coef_norm_Macaca_mulatta, y = Coef_norm_Rattus_norvegicus, 
                fill = celltype)) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 0, color = 'black') +
  geom_vline(xintercept = 0, color = 'black') +
  geom_point(pch = 21, size = 2) +
  geom_errorbar(aes(ymin=Coef_norm_Rattus_norvegicus - Coef_norm_se_Rattus_norvegicus,
                    ymax=Coef_norm_Rattus_norvegicus + Coef_norm_se_Rattus_norvegicus), 
                width=2, alpha = .6) +
  geom_errorbarh(aes(xmin=Coef_norm_Macaca_mulatta - Coef_norm_se_Macaca_mulatta,
                     xmax=Coef_norm_Macaca_mulatta + Coef_norm_se_Macaca_mulatta), 
                 height=2, alpha = .6) +
  scale_fill_manual(values = col_celltype) +
  facet_wrap( ~trait, scales = 'fixed') +
  xlab('Macaca mulatta') + ylab('Rattus norvegicus') + 
  theme_classic(base_size = 6) 
print(pp)
dev.off()



plot_fn = here(FIGDIR, 'plots','Rat_Transgen_NAc.vs.Stauffer_Monkey_striatum.LDSC.ppt.withoutOutliers.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn)

pp = ggplot(data = df_wide %>% filter(!grepl('AD|OCD|Doz', trait)),
            aes(x = Coef_norm_Macaca_mulatta, y = Coef_norm_Rattus_norvegicus, 
                fill = celltype)) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 0, color = 'black') +
  geom_vline(xintercept = 0, color = 'black') +
  geom_point(pch = 21, size = 2) +
  geom_errorbar(aes(ymin=Coef_norm_Rattus_norvegicus - Coef_norm_se_Rattus_norvegicus,
                    ymax=Coef_norm_Rattus_norvegicus + Coef_norm_se_Rattus_norvegicus), 
                width=2, alpha = .6) +
  geom_errorbarh(aes(xmin=Coef_norm_Macaca_mulatta - Coef_norm_se_Macaca_mulatta,
                     xmax=Coef_norm_Macaca_mulatta + Coef_norm_se_Macaca_mulatta), 
                 height=2, alpha = .6) +
  scale_fill_manual(values = col_celltype) +
  facet_wrap( ~trait, scales = 'fixed') +
  xlab('Macaca mulatta') + ylab('Rattus norvegicus') + 
  theme_classic(base_size = 6) 
print(pp)
dev.off()









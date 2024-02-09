library(ArchR)
library(tidyverse)
library(swfdr)
library(data.table)
library(writexl)
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

##############################################################################
## 1) merge the differential motif deviations across celltypes and groups
file_path <- here("data/differential_Motifs_data")
rds_files <- list.files(file_path, pattern = "differential", full.names = TRUE)
names(rds_files) = basename(rds_files) %>% 
  str_replace('differential_MotifMatrix\\.', '') %>% 
  str_replace('\\.rds', '') 

df = rds_files %>% lapply(function(file){
  # read in list of SummarisedExperiments
  df = readRDS(file) %>% 
    lapply(function(se){
      # extract the differential values from SE
      df = assays(se) %>% as.list() %>% lapply(as.data.frame) %>% bind_cols()
      names(df) = names(assays(se))
      # combine w/ motif data
      df = cbind(rowData(se) %>% as.data.frame(), df) %>% 
        mutate(TF = ss(name, '_')) %>% relocate(TF, .after = 'name') %>% 
        dplyr::select(-c(name, idx, seqnames))
      return(df)
    }) %>% rbindlist(idcol = 'group')
  return(df)
}) %>% rbindlist(idcol = 'celltype')

## recalculate FDR across all comparisons
df = df %>% mutate(
  group = group %>% 
    str_replace_all('Cocaine_v_Saline', 'Coc_vs_Sal') %>% 
    str_replace_all('Meth_v_Cocaine', 'Met_vs_Coc') %>% 
    str_replace_all('Meth_v_Saline', 'Met_vs_Sal'),
  FDR =  lm_qvalue(Pval, X=MeanBGD)$q
)
df %>% filter(FDR < 0.05) %>% count(celltype, group) %>% arrange(desc(n))

out_fn = here("data/differential_Motifs_data", 'combined_MotifMatrix.rds')
saveRDS(df, out_fn)

##############################################################################
## 2) merge the differential motif deviations across celltypes and groups
corMat = here('data/tidy_data/differential_expression/rdas/correlationMatrix.rds') %>% 
  readRDS() %>% as.data.frame() %>% dplyr::select(-MotifMatrix_name) %>% 
  dplyr::rename('TF' = 'GeneExpressionMatrix_name') %>% 
  mutate(isVariable = maxDelta > quantile(maxDelta, .5),
    TFRegulator = case_when( cor > 0.25 & padj < 0.01 & isVariable ~'Yes', T ~ "No")) %>% 
  dplyr::select(-c(GeneExpressionMatrix_seqnames:MotifMatrix_matchName))

corMat %>% count(TFRegulator)

corMat_signif = corMat %>% filter(TFRegulator == 'Yes')

## filter TFs that are correlated w/ expression
## recalculate the FDR of deviations across cell types 
df_cor = df %>% filter(TF %in% corMat_signif$TF) %>% 
  mutate(FDR =  lm_qvalue(Pval, X=MeanBGD)$q)
df_cor %>% filter(FDR < 0.05) %>% count(celltype, group) %>% arrange(desc(n))

out_fn2 = here("data/tidy_data", 'differential_expression', 'rdas', 'combined_MotifMatrix_TFRegulator.rds')
saveRDS(df_cor, out_fn2)


xcl_fn2 = here("data/tidy_data", 'differential_expression', 'tables', 'differential_correlatedMotifs_3-way_sires.xlsx')
write_xlsx(df_cor, xcl_fn2)

##############################################################################
## 3) compare the differential expression across groups
df_wide = df_cor %>%
  pivot_wider(id_cols = c(celltype, TF), names_from = 'group', 
              values_from = c(Mean:MeanBGD), values_fn = mean) %>% 
  mutate(
    ## calculate strength of DEG of stimulant vs. 
    effect_Met = -log10(Pval_Met_vs_Sal) * sign(MeanDiff_Met_vs_Sal),
    effect_Coc = -log10(Pval_Coc_vs_Sal) * sign(MeanDiff_Coc_vs_Sal),
    ## calculate an interaction score, weighted by MeanDifference
    interaction_MvC = (effect_Met - effect_Coc), 
    ## see if the differences between effects are concordant w/ MvC condition
    isConcordant = sign(interaction_MvC) == sign(MeanDiff_Met_vs_Coc)) %>% 
  arrange(desc(abs(interaction_MvC)))

out_fn3 = here("data/differential_Motifs_data", 'combined_MotifMatrix_TFRegulator_interaction.rds')
saveRDS(df_wide, out_fn3)

xcl_fn3 = here("data/tidy_data", 'differential_expression', 'tables', 'differential_correlatedMotifs_interactions.xlsx')
write_xlsx(df_wide, xcl_fn3)

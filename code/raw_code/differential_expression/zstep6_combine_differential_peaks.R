library(ArchR)
library(tidyverse)
library(swfdr)
library(data.table)
library(writexl)
library(qvalue)
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

##############################################################################
## 1) merge the differential Peak deviations across celltypes and groups
file_path <- here("data/differential_Peaks_data")
rds_files <- list.files(file_path, pattern = "differential", full.names = TRUE)
names(rds_files) = basename(rds_files) %>% 
  str_replace('differential_PeakMatrix\\.', '') %>% 
  str_replace('\\.rds', '') 

df = rds_files %>% lapply(function(file){
  # read in list of SummarisedExperiments
  df = readRDS(file) %>% 
    lapply(function(se){
      # extract the differential values from SE
      df =assays(se) %>% as.list() %>% lapply(as.data.frame) %>% bind_cols()
      names(df) = names(assays(se))
      # combine w/ Peak data
      df = cbind(rowData(se) %>% as.data.frame(), df) %>% 
        mutate(name = paste0(seqnames, ':',start, '-', end )) %>% 
        relocate(name, .before = 'Log2FC') %>% 
        dplyr::select(-c(idx))
      return(df)
    }) %>% rbindlist(idcol = 'group')
  return(df)
}) %>% rbindlist(idcol = 'celltype')

## recalculate FDR across all comparisons
## exclude peaks where no difference b/t groups b/c no data
df = df %>% filter(MeanBGD > 0) %>% mutate(
  group = group %>% 
    str_replace_all('Cocaine_v_Saline', 'Coc_vs_Sal') %>% 
    str_replace_all('Meth_v_Cocaine', 'Met_vs_Coc') %>% 
    str_replace_all('Meth_v_Saline', 'Met_vs_Sal'),
  FDR = qvalue(Pval)$q
) %>% arrange(FDR, Pval)
df %>% filter(FDR < 0.05) %>% count(celltype, group) %>% arrange(desc(n))

out_fn = here("data/differential_Peaks_data", 'combined_PeakMatrix.rds')
saveRDS(df, out_fn)

##############################################################################
## 2) extract peaks that correlate with gene expression to subset diff peaks
proj = here('data/tidy_data/ArchRProjects/Rat_Transgen_NAc_multiome') %>% 
  loadArchRProject(showLogo = F)
p2g <- getPeak2GeneLinks(proj, corCutOff = 0.25, resolution = 1, returnLoops = FALSE)
peaks = metadata(p2g)$peakSet %>% as.data.frame() %>% 
  mutate(name = paste0(seqnames, ':',start, '-', end ))
genes = metadata(p2g)$geneSet %>% mcols() %>% as.data.frame() %>% dplyr::select(-idx)
df_p2g = cbind( peaks[p2g$idxATAC,], 'gene' = genes[p2g$idxRNA,], as.data.frame(p2g)) %>% 
  dplyr::select(-starts_with('idx'))

df_p2g2 = df_p2g %>% arrange(FDR) %>% group_by(name) %>% 
  mutate(gene = paste(gene, collapse = ', ')) %>% ungroup() %>% 
  dplyr::select(-c(Correlation:VarQRNA))

### check if peaks same as diff peaks
table(peaks$name %in% df$name) # yes

## filter peaks that are correlated w/ expression
## recalculate the FDR of deviations across cell types 
df_cor = df %>% filter(name %in% df_p2g$name) %>% filter(MeanBGD > 0) %>% 
  mutate(FDR = qvalue(Pval)$q)
df_cor = cbind(df_cor, df_p2g2[match(df_cor$name, df_p2g2$name),'gene'])
df_cor %>% filter(FDR < 0.05) %>% count(celltype, group) %>% arrange(desc(n))

out_fn2 = here("data/tidy_data", 'differential_expression', 'rdas', 'combined_PeakMatrix_correlatedPeaks.rds')
saveRDS(df_cor, out_fn2)


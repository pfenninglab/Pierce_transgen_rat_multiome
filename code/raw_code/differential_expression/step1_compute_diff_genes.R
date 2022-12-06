## packages for data table processing 
library(here)
library(tidyverse)
library(data.table)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)
library(harmony)

## statistical genomics packages
library(countsplit) ## to control Type 1 errors
library(swfdr) ## to increase power of discovery

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

#######################################################
# 0) Seurat uses the future package for parallelization
plan("multicore", workers = 28)
options(future.globals.maxSize = 100 * 1024^3)
options(future.rng.onMisuse = 'ignore')

########################################
# 1) load in the full snRNA-seq dataset 
obj = here('data/tidy_data/Seurat_projects',
           "Rat_transgen_multiomeRNA_refined_all_SeuratObj_N5.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA') ## only the RNA data

##################################################################################
## 2) use countsplitting to control inflated p-values, http://arxiv.org/abs/2207.00554
## following https://anna-neufeld.github.io/countsplit/articles/seurat_tutorial.html
## useful pictoral description https://twitter.com/AnnaCNeufeld/status/1546598129637634048

set.seed(20221129)
## this subsets into 2 data, one set for re-clustering, one set for countspliting DEGs
## every cell is split into 2 cells, one in the countsplit data and one the train data
split <- countsplit(GetAssayData(obj, "counts"), epsilon=0.5)
Xtrain <- split$train %>% as("sparseMatrix")
Xcountsplit <- split$countsplit %>% as("sparseMatrix")
rm(split); gc()

all(colnames(Xcountsplit) == colnames(Xtrain))

## make the subsetted data for clustering
obj.train <- CreateSeuratObject(counts = Xtrain, min.cells = 3, min.features = 200, 
                                meta.data = obj[[]][colnames(Xtrain),])
Xtestsubset <- Xcountsplit[rownames(obj.train),colnames(obj.train)] ## create subset for later
dim(Xtestsubset)

## take the train subset of the data through the normalization and scaling 
obj.train <- obj.train %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

all.genes.train <- rownames(obj.train)

obj.train <- obj.train %>% 
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.ribo.mito") %>% 
  ## use Harmony to remove batch effects across samples
  RunPCA(verbose = TRUE) %>% 
  RunHarmony(group.by.vars = 'Sample', assay.use = "SCT", verbose = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  ## create data-based clustering
  FindNeighbors(dims = 1:30, verbose = TRUE) %>% 
  FindClusters(resolution = 2, algorithm = 2, verbose = TRUE) ## overcluster


## map the training clusterings w/ the annotation clusters
tt = with(obj.train[[]], table(seurat_clusters, cluster_rat))
relabel = setNames(colnames(tt)[apply(tt,1, which.max)], rownames(tt))

## just so the clusters were obtained mostly through clustering Xtrain
obj.train$cluster_rat2 = relabel[obj.train$seurat_clusters]

## make sure all re-clusterings have labels that are accounted for
length(unique(relabel)) == length(unique(obj.train$cluster_rat))

################################################################################################
## 3) add the countsplit data, to compute the DEGs according to cluster groupings of the training data
obj.train[['countsplit']] <- CreateAssayObject(counts=Xtestsubset)
obj.train <- obj.train %>%
  ## Normalize data w/ log
  NormalizeData(assay = 'countsplit') %>% ScaleData(assay = 'countsplit') %>% 
  ## use SCTransform to normalize data and regress out QC variables
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.ribo.mito", 
              assay = "countsplit", new.assay.name = 'countsplit_SCT') 

## save this just in case need to re-access the count-split object
obj.train %>% SaveH5Seurat(
  here('data/tidy_data/Seurat_projects',
       "Rat_transgen_multiomeRNA_refined_all_SeuratObj_N5_countsplit.h5Seurat"), 
  overwrite = T)

Idents(obj.train) = 'Sire'
df_cocVsal = FindConservedMarkers( 
  obj.train, ident.1 = 'Cocaine', ident.2 = 'Saline', 
  grouping.var = 'cluster_rat2', assay = "countsplit", slot = "data", 
  min.cells.group = 20, test.use = "wilcox", densify = T,
  min.pct=0, logfc.threshold = 0, verbose = TRUE)

df_methVsal = FindConservedMarkers( 
  obj.train, ident.1 = 'Methamphetamine', ident.2 = 'Saline', 
  grouping.var = 'cluster_rat2', assay = "countsplit", slot = "data", 
  min.cells.group = 20, test.use = "wilcox", densify = T,
  min.pct=0, logfc.threshold = 0, verbose = TRUE)

df_methVcoc = FindConservedMarkers( 
  obj.train, ident.1 = 'Methamphetamine', ident.2 = 'Cocaine', 
  grouping.var = 'cluster_rat2', assay = "countsplit", slot = "data", 
  min.cells.group = 20, test.use = "wilcox", densify = T,
  min.pct=0, logfc.threshold = 0, verbose = TRUE)

dir.create(here('data/tidy_data/differential_expression/rdas'), showWarnings = F)
save(df_cocVsal, df_methVsal, df_methVcoc, 
     file = here('data/tidy_data/differential_expression/rdas/diff_gene_3-way_countsplit.rda'))

#######################################################################################
# 4) combine the dataframes of differential expression by wilcoxon rank sum countsplit together

## helper variables for reshaping dataframe
celltypeMatch = obj.train$cluster_rat2 %>% unique() %>% paste0('^',.,'_') %>% paste(collapse = '|')
metricMatch = names(df_methVcoc) %>% grep(celltypeMatch, ., value = T)  %>%  
  str_replace(celltypeMatch, '') %>% unique() %>% paste0('_', ., '$') %>% paste(collapse = '|')

## gene average expression of each gene in each cell type, across Sire conditions
avg_expr = AverageExpression(obj.train, assay = "countsplit_SCT", group.by = "cluster_rat2")[[1]]
avg_expr = avg_expr %>% as.data.frame() %>% rownames_to_column('gene') %>% 
  pivot_longer(-gene, names_to = 'celltype', values_to = 'AveExpr')

## reshape the DEG table
df = list( 'Coc_vs_Sal' = df_cocVsal,
           'Met_vs_Sal' = df_methVsal,
           'Met_vs_Coc' = df_methVcoc) %>% 
  lapply(rownames_to_column, 'gene') %>% rbindlist(idcol = 'group') %>% 
  pivot_longer(cols = -c(group, gene)) %>% 
  mutate(metric = str_replace_all(name, celltypeMatch, ''),
         celltype = str_replace_all(name, metricMatch, '')) %>% 
  dplyr::select(-name) %>% filter(!metric %in% c('max_pval', 'minimump_p_val')) %>% 
  ## combine the percent expressed in XX group number
  pivot_wider(values_from = 'value', names_from = 'metric')

## use average expression matrix to model power in DEG analyses
df = df %>% inner_join(avg_expr) %>% 
  ## use SWFDR to increase power of detecting DEGs based on avg expression covariate
  ## https://pubmed.ncbi.nlm.nih.gov/30581661/
  ## correction is done across all pairwise grouping and cell types
  mutate(p_val_adj = lm_qvalue(p_val, X=AveExpr)$q, 
         p_val_bonf = p.adjust(p_val, 'bonferroni'))

################################################################
# 5) nominate DEGs by the degree these cells are expressed
df_wide = df %>% dplyr::select(-c(pct.1, pct.2)) %>% 
  relocate(c(gene, celltype, AveExpr), .before = everything()) %>% 
  pivot_wider(id_cols = gene:AveExpr, names_from = 'group', 
              values_from = c(p_val:p_val_bonf)) %>% 
  mutate(
    ## calculate strength of DEG of stimulant vs. 
    effect_Met = -log10(p_val_Met_vs_Sal) * avg_log2FC_Met_vs_Sal,
    effect_Coc = -log10(p_val_Coc_vs_Sal) * avg_log2FC_Coc_vs_Sal,
    ## calculate an interaction score, weighted by AveExpr
    interaction_MvC = (effect_Met - effect_Coc)/ abs(effect_Met + effect_Coc)) %>% 
  dplyr::select(-c(effect_Met, effect_Coc)) %>% 
  arrange(desc(abs(interaction_MvC)))

## list of dataframes of DEGs by cell type
df_list = df_wide %>% split(f = .$celltype)
names(df_list) = names(df_list) %>% make.names()

## candidate list of DEGs that are different b/c Meth vs. Sal and Coc vs. Sal
alpha = 0.05
df_interact = df_wide %>% 
  ## one of the vs. saline condition must be significant
  filter(p_val_bonf_Coc_vs_Sal < alpha | p_val_bonf_Met_vs_Sal < alpha |
           p_val_bonf_Coc_vs_Sal < alpha) %>% 
  ## the interaction meth vs. coc must be significant
  filter(abs(interaction_MvC)>1)
  

######################################################
# 6) look at overall trends and export the DEG tables
## do some counting of DEGs
df %>% group_by(group, celltype) %>% 
  summarise(numDEG = sum(p_val_adj < alpha & abs(avg_log2FC) > .2) ) %>% 
  pivot_wider(names_from = 'group', values_from = 'numDEG')

df %>% group_by(group, celltype) %>% 
  summarise(numDEG = sum(p_val_bonf < alpha & abs(avg_log2FC) > .2) ) %>% 
  pivot_wider(names_from = 'group', values_from = 'numDEG')

table(df_interact$celltype)
df_interact %>% pull(gene) %>% unique() %>% sort() %>% paste(collapse = ', ')

saveRDS(df_wide, file = here('data/tidy_data/differential_expression/rdas/diff_gene_3-way_countsplit.rds'))
dir.create(here('data/tidy_data/differential_expression/tables'), showWarnings = F)
df_list %>% writexl::write_xlsx(here('data/tidy_data/differential_expression/tables',
                                     'differential_expression_3-way_sires.xlsx'))
df_interact %>% writexl::write_xlsx(here('data/tidy_data/differential_expression/tables',
                                         'differential_expression_top_interaction_MvC.xlsx'))


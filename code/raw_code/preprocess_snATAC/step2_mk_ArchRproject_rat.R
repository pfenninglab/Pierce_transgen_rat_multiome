suppressMessages(library(ArchR))
library(here)
library(tidyverse)
library(parallel)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = min(length(fragments_fn), round(parallel::detectCores()/4)))
GENOMEDIR='/home/bnphan/resources/genomes/rn7'
load(file.path(GENOMEDIR,'rn7_liftoff_mm10NcbiRefSeq_ArchR_annotations.rda'))

#############################################################################
### make that arrow file from bgzipped fragments file (tsv.gz or bed.gz) ####
ArrowFiles <- here('data/raw_data/arrow') %>%
  list.files(pattern = '.arrow$', full.names = T)

## drop sample 6 b/c of low QC
ArrowFiles = ArrowFiles[!grepl('6C', ArrowFiles)]

###############################
### make the ArchR Project ####
dir.create(here('data/tidy_data/ArchRProjects'), showWarnings = F, recursive = T)
proj = ArchRProject( 
  ArrowFiles = ArrowFiles, copyArrows = TRUE,
  geneAnnotation = geneAnnotation, #must be the custom rn7 version
  genomeAnnotation = genomeAnnotation, #must be the custom rn7 version
  outputDirectory = here("data/tidy_data/ArchRProjects/Rat_Transgen_NAc_scATAC"))

## filter doublets
proj = filterDoublets( proj, cutEnrich = .5, cutScore = -Inf, filterRatio = 1)

## increase unique fragments cutoff to 10^3.5, remove cluster of low QC cell 
idxSample <- BiocGenerics::which(proj$nFrags > 10^3.6)
cellsSample <- proj$cellNames[idxSample]
proj = subsetCells(ArchRProj = proj, cellNames = cellsSample)

## Re-create the metadata from sample names
pd = data.frame(cellNames = getCellNames(proj)) %>% 
  mutate(sampleName = ss(cellNames, '#'), 
         Sire = case_when(grepl('C', sampleName) ~ 'Cocaine', 
                          grepl('M', sampleName) ~ 'Methamphetamine', 
                          TRUE ~ 'Saline'),
         ) %>% dplyr::select(-c(sampleName)) %>% 
  column_to_rownames('cellNames')

for(col in names(pd)){
  proj = addCellColData(
    ArchRProj = proj,
    data = pd[proj$cellNames,col],
    name = col,
    cells = getCellNames(proj),
    force = TRUE
  )
}
table(proj$Sample, proj$Sire)

## save project
proj = saveArchRProject(proj)


suppressMessages(library(ArchR)) ## >= v1.0.2
library(here)
library(tidyverse)
library(parallel)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

###############################################
### get the fragments file from snATAC-seq ####
fragments_fn = here('data/raw_data/CellRanger_outs') %>%
  list.files(pattern = 'atac_fragments.tsv.gz$', full.names = T, recursive = T)
names(fragments_fn) = fragments_fn %>% ss('/', 9)

########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = min(length(fragments_fn), round(parallel::detectCores()/4)))
GENOMEDIR='/home/bnphan/resources/genomes/rn7'
# contains `geneAnnotation` and `genomeAnnotation` objects
load(file.path(GENOMEDIR,'rn7_liftoff_mm10NcbiRefSeq_ArchR_annotations.rda'))

#############################################################################
### make that arrow file from bgzipped fragments file (tsv.gz or bed.gz) ####
ArrowFiles <- createArrowFiles( inputFiles = fragments_fn, 
                                sampleNames = names(fragments_fn), 
                                minTSS = 2, #Dont set this too high because you can always increase later
                                minFrags = 1000, addTileMat = TRUE,
                                addGeneScoreMat = TRUE,
                                geneAnnotation = geneAnnotation, #must be the custom rn7 version
                                genomeAnnotation = genomeAnnotation, #must be the custom rn7 version
                                force = TRUE)

doubleScores <- addDoubletScores(input = ArrowFiles, 
                                 k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                                 LSIMethod = 1,  knnMethod = "UMAP") #Refers to the embedding to use for nearest neighbor search with doublet projection.

############################################
### move Arrow files and QC to data dir ####
rsync = 'rsync -Paq --remove-source-files'
files = c(paste0(names(fragments_fn),'.arrow'), 'QualityControl', 'ArchRLogs')
dir = here("data/raw_data/arrow")
thecall = paste(rsync, files, dir)
sapply(thecall, system)


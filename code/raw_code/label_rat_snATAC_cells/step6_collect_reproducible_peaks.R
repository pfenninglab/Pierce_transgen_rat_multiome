ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer))
library(tidyverse)
library(here)

source(here('code/final_code/hal_scripts/narrowPeakFunctions.R'))

##############################
### read in ArchR project ####
DATADIR='data/tidy_data'
CODEDIR='code/raw_code/label_rat_snATAC_cells'

##############################
### read in ArchR project ####
DATADIR='data/tidy_data'
LABEL='Rat_Transgen_NAc'; GENOME = 'rn7'; 

#######################################################
# get the reproducible peaks across Cluster2 cell types
peak_rds_fn = c(list.files(path = here(DATADIR,'ArchRProjects', 
                                       'Rat_Transgen_NAc_scATAC_clusterRat','PeakCalls'), 
                         full.names = T, pattern = '.rds'), 
                list.files(path = here(DATADIR,'ArchRProjects', 
                                       'Rat_Transgen_NAc_scATAC_clusterMacaque','PeakCalls'), 
                           full.names = T, pattern = '.rds'))
names(peak_rds_fn) = ss(basename(peak_rds_fn),'-reproduciblePeaks.gr.rds')
peakList = lapply(peak_rds_fn, readRDS)

# label the summit and peak name using rn7 coordinates
peakList = lapply(peakList, addSummitCenter)
peakList = lapply(peakList, nameNarrowPeakRanges, genome = GENOME)
peakList = lapply(peakList, sort)

###############################################
# create directory and narrowPeak file names ##
PEAKDIR2=here('data/raw_data','peak_rn7')
system(paste('mkdir -p',PEAKDIR2))
narrowPeak_rn7_fn = here(PEAKDIR2, paste0(LABEL, '.', names(peakList), 
                                             '.rn7.narrowPeak.gz'))
# write peaks to narrowPeak file
outRanges = mapply(write_GRangesToNarrowPeak,
                   gr = peakList, file = narrowPeak_rn7_fn, genome = GENOME)

############################################################
## liftover peaks to rn6 (the version in Cactus hal file) ##
chainFile =here("/home/bnphan/resources/liftOver_chainz", 'rn7ToRn6.over.chain')
peakList_rn6 = lapply(peakList, liftOver_narrowPeak, chainFile = chainFile)

PEAKDIR=here('data/raw_data','peak_cluster_macaque')
dir.create(PEAKDIR, showWarnings = F, recursive = T)
narrowPeak_fn = here(PEAKDIR, paste0(LABEL, '.', names(peakList_rn6), 
                                     '.UCSCRn6.narrowPeak.gz'))
# write peaks to narrowPeak file
outRanges = mapply(write_GRangesToNarrowPeak,
                   gr = peakList_rn6, 
                   file = narrowPeak_fn, genome = GENOME)

sapply(outRanges, function(gr){
  tmp = gr %>% as.data.frame() %>%
    mutate(col = paste0(seqnames, start, end, peak, name)) %>% 
    filter(duplicated(col))
  return(nrow(tmp))
})



# take care of w/ Carson's 

#####################################
# halLiftOver and HALPER the peaks ##
# map to human
halmapper_script = '/home/bnphan/src/atac_data_pipeline/scripts/halper_map_peak_orthologs.sh'
system('mkdir -p logs')
sbatch = 'sbatch -p pfen1 -w compute-1-11 --mem 10G'
SOURCE_SPECIES = 'Rattus_norvegicus'
TARGET_SPECIES = c('Homo_sapiens,Mus_musculus,Macaca_mulatta')
target_species = paste('-t', TARGET_SPECIES)
source_species = paste('-s', SOURCE_SPECIES)
outdir = paste('-o', here('data/raw_data', 'halper'))
peak_files = paste('-b',narrowPeak_fn)
halLiftover = paste('--halPath', '/projects/pfenninggroup/install/halLiftover')

# paste the parameter calls together
thecall1 = paste(sbatch, halmapper_script, source_species,
                rep(target_species, each= length(peak_files)),
                outdir, rep(peak_files, times = length(target_species)),
                halLiftover)

cat(thecall1, file= here(CODEDIR,'step7_map_to_other_species.sh'), sep = '\n')
system(paste('chmod u+x', here(CODEDIR, 'step7_map_to_other_species.sh')))
system('./step7_map_to_other_species.sh')



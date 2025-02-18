library(ArchR)
library(ArchRtoSignac)
library(SeuratDisk)
library(Seurat)
library(Signac)
library(here)
library(tidyverse)
library(rtracklayer)
library(future)
library(harmony)
library(EnsDb.Rnorvegicus.v79)
library(BSgenome.Rnorvegicus.UCSC.rn7)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

PLOTDIR='figures/exploratory/multiome_integration'
sapply(here(PLOTDIR, c('plots', 'tables', 'rdas')), dir.create, showWarnings = F, recursive = T)

#########################################################
# 0) Seurat uses the future package for parallelization
plan("multicore", workers = 12)
options(future.rng.onMisuse = 'ignore')
options(future.globals.maxSize = 60 * 1024^3)

celltypes =  c("D1-MSN", "D1/D3-MSN", "D1/D2-Hybrid-MSN", 
               "D1-NUDAP-MSN", "D1-ICj-MSN", "D2-MSN", "Pvalb-INT", "Sst-INT" ,
               "Astrocytes", "Microglia", "Oligos", "Oligos_Pre")

cols = paletteDiscrete(values = celltypes)


##############################
## 1) load in the Seurat object
obj = LoadH5Seurat(here::here('data/tidy_data/Seurat_projects', 
                              "Rat_transgen_RNA-ATAC_SeuratObj_N5.h5Seurat"))

Idents(obj) <- "cluster_rat"
levels(obj) <- celltypes

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Rnorvegicus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "rn7"

# add the gene information to the object
Annotation(obj) <- annotations

# compute gene activities
gene.activities <- GeneActivity(brain)

obj[['RNA']] <- CreateAssayObject(counts = gene.activities)
obj <- NormalizeData(
  object = obj,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(obj$nCount_RNA)
)


#################################
## 2) compute peak2gene correlation

markerGenes = c(
  'RBFOX3', # NeuN-neurons'
  'GAD2', # Gaba-ergic neurons
  'PPP1R1B', # DARPP-32 MSNs
  'DRD1', # direct MSNs
  'OPRM1', # NUDAP
  'RXFP1', # D1/D2 hybrid
  'DRD3', # ICj 
  'DRD2',  # indirect MSNs
  'LHX6', # MGE interneurons
  'AQP4', # Astrocytes
  'CX3CR1', # Microglia
  'MOG', # Oligodendrocyte
  'PDGFRA' # OPC
) %>% str_to_title()

obj <- LinkPeaks( obj, peak.assay = "ATAC", expression.assay = "MAGIC_RNA",
                   genes.use = c('Camk4', 'Drd1', 'Drd2'))

Signac:::FindRegion(obj, 'Drd1')

######################################################
## 3) construct the multipanel multiomic track plots
idents.plot <- celltypes %>% str_subset('MSN')

obj2 = subset(obj, cluster_rat %in%idents.plot & cluster_rat != 'D1-ICj-MSN' )

for (gene in  c('Camk4', 'Drd1', 'Drd2')){
  p1 <- CoveragePlot( 
    obj2, group_by = 'cluster_rat', 
    region = gene, features = gene, expression.assay = "MAGIC_RNA",
    annotation = T, scale.factor = 1e7, ymax = 15,
    extend.upstream = 10000, extend.downstream = 10000) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8) # Rotate x-axis labels
    )
  
  plot_fn = here(PLOTDIR, 'plots', paste0('fixed_trackplot_', gene, '.neurons.pdf'))
  pdf(plot_fn, width = 5, height =3)
  print(p1 & scale_fill_manual(values = cols))
  dev.off()
}

# in2mm = 25.4
# fig2_diagonal_matrix_dotplot_fn = 
#   here(PLOTDIR, 'plots', 'diagonal_matrix_dotplot.all.pdf')


# # https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf
# pdf(fig2_diagonal_matrix_dotplot_fn, width = 110/in2mm, height =  70/in2mm)
# DotPlot( obj, features = markerGenes, cols =c('lightgray', 'navyblue'),
#          cluster.idents = F, scale = T, scale.by = "radius", assay = 'RNA') +
#   theme_bw() + scale_y_discrete(limits = rev) +
#   scale_x_discrete(position = "top") +
#   guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5),
#          shape = guide_legend(title.position="top", title.hjust = 0.5))+
#   theme(legend.position = 'bottom', plot.title= element_blank(),
#         axis.text.x = element_text(angle = -30, vjust = 0, hjust=1),
#         axis.title = element_blank(), 
#         legend.text = element_text(size = 6), 
#         legend.title = element_text(size = 6),
#         legend.key.size = unit(1, "mm"),
#         legend.spacing.x = unit(1, 'mm'))
# dev.off()



# Load required libraries
library(Seurat)        # Assuming you are working with Seurat
library(rtracklayer)   # For liftOver functionality
library(GenomicRanges) # For handling genomic data

# Assuming 'seurat_obj' is your Seurat object
# Extract genomic coordinates from rownames
genomic_coords <- rownames(obj)

# Take the first 100 genomic coordinates
genomic_coords_subset <- head(genomic_coords, 100)

# Create a GRanges object for rn6 coordinates
rn6_regions <- GRanges()

for (coord in genomic_coords_subset) {
  # Split the coordinates
  parts <- unlist(strsplit(coord, "-"))
  
  # Check if the format is correct (i.e., "chrN-start-end")
  if (length(parts) == 3) {
    chr <- parts[1]
    start <- as.numeric(parts[2])
    end <- as.numeric(parts[3])
    
    # Create GRanges object for the rn6 region
    region <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
    rn6_regions <- c(rn6_regions, region)
  }
}

# Load the liftOver chain file
chain_file <- "/home/csriniv1/resources/rn6ToRn7.over.chain.gz"  # Adjust the path to your downloaded chain file
chain <- import.chain(chain_file)

# Perform the liftover
lifted_regions <- liftOver(rn6_regions, chain)

# Check if any coordinates were successfully lifted
if (length(lifted_regions) > 0) {
  # Create a data frame to compare original and lifted coordinates
  comparison <- data.frame(
    Original = as.character(rn6_regions),
    Lifted = as.character(unlist(lifted_regions))
  )
  
  # Check if coordinates have changed
  comparison$Changed <- comparison$Original != comparison$Lifted
  
  # Print the comparison
  print(comparison)
  
  # Summary
  if (any(comparison$Changed)) {
    cat("Some coordinates have changed during liftover to rn7.\n")
  } else {
    cat("No coordinates changed during liftover to rn7.\n")
  }
} else {
  cat("No coordinates were successfully lifted to rn7.\n")
}
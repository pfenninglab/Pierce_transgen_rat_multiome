# Pierce_transgen_rat_multiome
Pierce lab rat nucleus accumbens single cell multiome analyses in transgenerational model of substance use predisposition

# Outline of of the main code:
These code folders are organized in roughly order of operations. The snRNA-seq data are clustered first, then they are used to aid snATAC-seq clustering. The differentially expressed gene analyses are done after the snRNA-seq clustering. The LD score regression (GWAS enrichment) analyses are done afte rthe snATAC-seq clustering.  
[code/raw_code/preprocess_snATAC](code/raw_code/preprocess_snATAC)
- code to make snATAC and snRNA raw files

[code/raw_code/label_rat_snRNA_cells](code/raw_code/label_rat_snRNA_cells)
- Seurat v4 clustering of the data an guided annotations with a reference dataset
- produces the `refined` clustering annotations based on this dataset `cluster_rat` or the reference He, Kleyman et al. dataset `cluster_macaque`.

[code/raw_code/label_rat_snATAC_cells](code/raw_code/label_rat_snATAC_cells)
- ArchR v1.0.2 clustering of the data borrowing the same cell barcodes of the snRNA-seq clustered data
- produces the cleaned up ArchR projects with the `cluster_rat` or `cluster_macaque` labels propagated across snATAC cells with or without an snRNA expression profile.
- also produces the reproducible peaks (peak found in N>=2 samples) in rn7 and rn6 coordinates
- maps these peaks to relevant species through a [240 mammals CACTUS MSA alignment v2](https://cglgenomics.ucsc.edu/data/cactus/) (hg38, mm10, rheMac8)

[code/raw_code/differential_expression](code/raw_code/differential_expression)
- differential expression using countsplitting
- produces the DEG lists for `Cocaine vs. Saline`, `Methamphetamine vs. Saline` and `Methamphetamine vs. Cocaine` pairwise comparisons.

[code/raw_code/ldsc_rat_snATAC](code/raw_code/ldsc_rat_snATAC)
- Stratified Linkage-Disequilibrium Score Regression (S-LDSC) analyses using the peaks mapped to hg38 coordinates


# Outline of of the supporting code directories:
[code/final_code/HeKleyman2021_macaque_striatum_data_processing](code/final_code/HeKleyman2021_macaque_striatum_data_processing)
- Downloads the He, Kleyman et al. snRNA-seq reference dataset. Subset the data to the NAc and recluster using sctransform and Seurat v4 integration methods
- Maps the monkey orthologs human genes used in this dataset to map to 1-1 mouse gene orthologs as defined by ENSEMBL.
- produces the files `GSE167920_Results_full_nuclei_processed_final_NAc.mmGenes.h5Seurat` and `GSE167920_Results_MSNs_processed_final_NAc.mmGenes.h5Seurat` used in the snRNA-seq guided cell type annotations

[code/final_code/hal_scripts/narrowPeakFunctions.R](code/final_code/hal_scripts/narrowPeakFunctions.R)
- Wrapper functions to facilitate mapping narrowPeak files between genome versions and cleaning up ArchR reproducible peaks

#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --exclude=compute-1-11,compute-1-12
##SBATCH --time 3-00:00:00
#SBATCH --job-name=grnboost
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=23G
#SBATCH --error=logs/pyscenic_grnboost_%A_%a.txt
#SBATCH --output=logs/pyscenic_grnboost_%A_%a.txt
#SBATCH --array=50-100

###################
## 1) set up paths 
PROJDIR=/projects/pfenninggroup/singleCell/Pierce_transgen_rat_multiome
DATADIR=$PROJDIR/data/tidy_data/scenic_grn
TF_LIST=/home/bnphan/resources/SCENIC/allTFs_mm.txt
DBNAME=/home/bnphan/resources/SCENIC/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
ANNOT=/home/bnphan/resources/SCENIC/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
mkdir -p $DATADIR/grnboost2
cd $PROJDIR/code/raw_code/scenic_grn

source ~/.bashrc; mamba activate pyscenic
rm -rf /tmp/*

########################################################################
## 1) run the GRNBoost2 w/ pysenic 100 times with the subset of CTL pts
N=`expr $SLURM_ARRAY_TASK_ID % 12 + 1`
ITER=`expr $SLURM_ARRAY_TASK_ID / 12 + 1`
EXPRMAT=$(ls $DATADIR/loom/* | sed "${N}q;d")
LABEL=$(basename $EXPRMAT .loom | sed 's/.*\.//g')

pyscenic ctx \
-o OUTPUT \
--mode custom_multiprocessing \
--annotations_fname ${ANNOT} \
--num_workers 8 \
--thresholds 0.90 \
--expression_mtx_fname ${EXPRMAT} \
--cell_id_attribute "CellID" \
--gene_attribute "Gene" \
module_fname ${DBNAME} 


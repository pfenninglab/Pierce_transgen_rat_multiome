#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --exclude=compute-1-11,compute-1-5,compute-1-35,compute-4-9
#SBATCH --time 3-00:00:00
#SBATCH --job-name=grnboost
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=23G
#SBATCH --error=logs/pyscenic_grnboost_%A_%a.txt
#SBATCH --output=logs/pyscenic_grnboost_%A_%a.txt
#SBATCH --array=25-99%5

###################
## 1) set up paths 
PROJDIR=/projects/pfenninggroup/singleCell/Pierce_transgen_rat_multiome
DATADIR=$PROJDIR/data/tidy_data/scenic_grn
TF_LIST=/home/bnphan/resources/SCENIC/allTFs_mm.txt
mkdir -p $DATADIR/grnboost2
cd $PROJDIR/code/raw_code/scenic_grn

source ~/.bashrc; mamba activate pyscenic; 

########################################################
## 1) run the GRNBoost2 w/ pysenic N times for this file
N=`expr $SLURM_ARRAY_TASK_ID % 5 + 1`
ITER=`expr $SLURM_ARRAY_TASK_ID / 5 + 1`
EXPRMAT=$(ls $DATADIR/loom/* | sed "${N}q;d")
LABEL=$(basename $EXPRMAT .loom | sed 's/.*\.//g')

OUTFILE=$DATADIR/grnboost2/Rat_transgen_multiomeRNA_GRNBoost2.${LABEL}.${ITER}.tsv
echo "Working on Sample ID:${LABEL} for the ${ITER}th run."
if [[ ! -f $OUTFILE || $(cat $OUTFILE) == '' ]]; then
arboreto_with_multiprocessing.py \
-o ${OUTFILE} -m grnboost2 \
--seed ${SLURM_ARRAY_TASK_ID} \
--num_workers 8 \
--cell_id_attribute "CellID" \
--gene_attribute "Gene" \
${EXPRMAT} ${TF_LIST}
fi

exit
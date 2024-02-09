#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=interactive
#SBATCH --time 8:00:00
#SBATCH --job-name=combine
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --dependency=afterok:3311939
#SBATCH --mem=10G
#SBATCH --error=logs/combine_grnboost_%A_%a.txt
#SBATCH --output=logs/combine_grnboost_%A_%a.txt
#SBATCH --array=1-5

###################
## 1) set up paths 
PROJDIR=/projects/pfenninggroup/singleCell/Pierce_transgen_rat_multiome
DATADIR=$PROJDIR/data/tidy_data/scenic_grn
mkdir -p $DATADIR/grnboost2
cd $PROJDIR/code/raw_code/scenic_grn

########################################################
## 1) find replicate runs of GRNboost for a given sample
EXPRMAT=$(ls $DATADIR/loom/* | sed "${SLURM_ARRAY_TASK_ID}q;d")
LABEL=$(basename $EXPRMAT .loom | sed 's/.*\.//g')
REP=$(ls $DATADIR/grnboost2/Rat_transgen_multiomeRNA_GRNBoost2.${LABEL}.*.tsv| wc -l|cut -f1)

##################################################
## merge the TF-target importances across samples
REP=12 ## number of samples
COUNT=$DATADIR/grnboost2/Combined_Rat_transgen_multiomeRNA_GRNBoost2.tsv
echo -e "TF\ttarget\timportance" > $COUNT
awk FNR!=1 $DATADIR/grnboost2/Avg_GRNBoost2.*.tsv | \
awk -F '\t' '
BEGIN { OFS = "\t"} 
NR>1{ $4 =$1"#"$2; arr[$4] += $3; count[$4] += 1} 
END{ print "TF\ttarget\timportance\n"}
{for (a in arr) { if (count[a] >= REP * .8){
print a OFS arr[a]/count[a]
}}}' | tr '#' '\t' >> $COUNT 

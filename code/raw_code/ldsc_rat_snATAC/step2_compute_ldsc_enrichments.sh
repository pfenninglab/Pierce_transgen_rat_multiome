#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time 1-0
#SBATCH --job-name=ldsc_rat
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/ldsc_rat_%A_%a.txt
#SBATCH --output=logs/ldsc_rat_%A_%a.txt
#SBATCH --array=1-74

# get the GWAS for this array job
SETDIR=/projects/pfenninggroup/singleCell/Pierce_transgen_rat_multiome
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
CODEDIR=${SETDIR}/code/raw_code/ldsc_rat_snATAC
DATADIR=${SETDIR}/data/tidy_data/ldsc_rat_snATAC
cd $CODEDIR; source activate ldsc

# get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas_list_sumstats.tsv )
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas_list_sumstats.tsv )
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas_list_sumstats.tsv )

OUTDIR=${DATADIR}/enrichments; mkdir -p $OUTDIR

#################################################################################
# run LD score regression over the Stauffer striatum cell type binary annotations
CTS_FN1=${DATADIR}/Rat_Transgen_NAc.${POP}_hg38.ldcts
if [[ ! -f "$OUTDIR/Rat_Transgen_NAc.${GWAS_Label}.${POP}.cell_type_results.txt" ]]; then
ldsc.py --ref-ld-chr-cts $CTS_FN1 \
--ref-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/baseline_v1.1/baseline_v1.1.${POP}. \
--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
--h2-cts $GWAS --out $OUTDIR/Rat_Transgen_NAc.${GWAS_Label}.${POP}
fi

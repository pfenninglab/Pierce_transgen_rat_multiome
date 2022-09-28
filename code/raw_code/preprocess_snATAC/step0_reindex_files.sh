PROJDIR=/projects/pfenninggroup/singleCell/Pierce_transgen_rat_multiome
cd $PROJDIR

for FILE in $PROJDIR/data/raw_data/CellRanger_outs/*/atac_fragments.tsv.gz; do
echo $FILE
tabix -p bed ${FILE}
done
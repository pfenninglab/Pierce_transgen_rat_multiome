PROJDIR=/projects/pfenninggroup/singleCell/Pierce_transgen_rat_multiome
cd $PROJDIR

for FILE in $PROJDIR/data/raw_data/CellRanger_outs/*/atac_fragments.tsv.gz; do
echo $FILE
tabix -p bed ${FILE}
done

## make a copy of the fragments to renamed w/ sample ID in file
for DIR2 in $PROJDIR/data/raw_data/fragments/*/outs; do
cd $DIR2
mv atac_fragments.tsv.gz fragments.tsv.gz
mv atac_fragments.tsv.gz.tbi fragments.tsv.gz.tbi
done


## make a copy of the fragments to renamed w/ sample ID in file
for DIR in $PROJDIR/data/raw_data/fragments/*/outs; do
DIR2=$(echo $DIR | sed 's/\/fragments\//\/CellRanger_outs\//g' |\
	sed 's/\/outs/\//g' )
rsync -Paq $DIR/*fragments* $DIR2
done


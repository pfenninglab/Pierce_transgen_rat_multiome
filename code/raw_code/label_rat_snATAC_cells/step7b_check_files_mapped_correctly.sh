###############################################
# 2) Remove early terminated halLiftover jobs #
function checkFile(){
HASFILE=TRUE
NUMLIFT=$(zcat $FILE | cut -f4 | cut -f1-2 -d':'| sort | uniq | wc -l | cut -d ' ' -f1)
gzip -t $FILE
if [[ $(echo $?) != '0' ]]; then 
# if gzipped file is corrupt
echo "gzipped file is corrupt: $(basename $FILE)."
HASFILE=FALSE; rm $FILE
elif [[ $NUMLIFT < 20 ]]; then 
echo "too few source chromosomes in $FILE."
HASFILE=FALSE; rm $FILE; fi
}

#### remove files that didn't have complete annotations
TOLOOK=$(ls *.HALPER.narrowPeak.gz)
for FILE in $TOLOOK; do
NAME=$(basename $FILE | sed 's/.HALPER.narrowPeak.gz//g' |sed 's/Homo_sapiensTo/mappableTo./g')
# check if mapping is right
checkFile
if [[ $HASFILE == "FALSE" ]]; then
echo "Removing $(basename $FILE .HALPER.narrowPeak.gz) files."
rm $(echo $FILE | sed 's/.HALPER.narrowPeak.gz/.halLiftover.sFile.bed.gz/g') $(echo $FILE | sed 's/.HALPER.narrowPeak.gz/.halLiftover.tFile.bed.gz/g')
fi
done
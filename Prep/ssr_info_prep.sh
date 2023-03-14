#!/bin/bash

# script to convert SSR info csv files to tab delimeted text files
# usage: bash ssr_info_prep.sh
# usage Ex: bash ssr_info_prep.sh

# retrieve ssr info path
infoInput=$(grep "info:" ../"InputData/inputPaths_ssr_pipeline.txt" | tr -d " " | sed "s/info://g")
# retrieve sequence info path
referenceInput=$(grep "reference:" ../"InputData/inputPaths_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")

# remove txt extension from input info file
infoOutput=$(echo $infoInput | sed 's/\.csv/\.txt/g')

# exclude the header of the ssr info file
# and exclude the lines beginning with empty cells that contain notes
# and fix any marker tags by removing excess > symbols
# and convert the delimeter from commas to tabs
# and retrieve the following fields:
# Marker name, start repeat, end repeat, Sequence showing primer sequences and repeat in BOLD
cat $infoInput | tail -n+2 | sed '/^,/d' | sed 's/>//g' | sed 's/,/\t/g' | cut -f1-3,9 > $infoOutput".tmp.txt"

# retrieve the marker ssr info
cat $infoOutput".tmp.txt" | cut -f1-3 > $infoOutput".tmpcol1.txt"

# retrieve the 62 marker sequences and exclude the final empty line
cat $infoOutput".tmp.txt" | cut -f4 | sed '$ s/$/'"$markerSeq"'/' | head -62 > $infoOutput".tmpcol2.txt"

# retrieve the missing marker sequence
markerSeq=$(cat $referenceInput | grep -A 1 "054-CM_015" | tail -1)

# add the missing marker sequence
echo $markerSeq >> $infoOutput".tmpcol2.txt"

# combine the marker ssr info and sequences columns
paste $infoOutput".tmpcol1.txt" $infoOutput".tmpcol2.txt" > $infoOutput

# clean up
rm $infoOutput".tmp.txt"
rm $infoOutput".tmpcol1.txt"
rm $infoOutput".tmpcol2.txt"

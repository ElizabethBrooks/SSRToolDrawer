#!/bin/bash

# script to convert SSR info csv files to tab delimeted text files
# usage: bash ssr_info_formatting.sh
# usage Ex: bash ssr_info_formatting.sh

# retrieve ssr info path
infoOutput=$(grep "ssrInfo:" ../"InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrInfo://g")
#infoOutput="/Users/bamflappy/GBCF/JRS/romero_Mar2023/Data/SSRinfo-2023.txt"
# retrieve sequence info path
referenceInput=$(grep "reference:" ../"InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")
#referenceInput="/Users/bamflappy/GBCF/JRS/romero_Mar2023/Data/finalset-2.fa"

# remove extension from input info file
infoInput=$(echo $infoOutput | sed 's/\.txt//g')

# remove hidden carriage return characters
# and exclude the header of the ssr info file
# and exclude the lines beginning with empty cells that contain notes
# and fix any marker tags by removing excess > symbols
# and convert the delimeter from commas to tabs
# and retrieve the following fields:
# Marker name, start repeat, end repeat, Sequence showing primer sequences and repeat in BOLD
cat $infoInput".csv" |  tr -d $'\r' | tail -n+2 | sed '/^,/d' | sed 's/>//g' | sed 's/,/\t/g' | cut -f1-3,9 > $infoInput".tmp.txt"

# retrieve the marker ssr info
cat $infoInput".tmp.txt" | cut -f1-3 > $infoInput".tmpcol1.txt"

# retrieve the 62 marker sequences and exclude the final empty line
cat $infoInput".tmp.txt" | cut -f4 | sed '$ s/$/'"$markerSeq"'/' | head -62 > $infoInput".tmpcol2.txt"

# retrieve the missing marker sequence
markerSeq=$(cat $referenceInput | grep -A 1 "054-CM_015" | tail -1)

# add the missing marker sequence
echo $markerSeq >> $infoInput".tmpcol2.txt"

# combine the marker ssr info and sequences columns
paste $infoInput".tmpcol1.txt" $infoInput".tmpcol2.txt" > $infoOutput

# clean up
rm $infoInput".tmp.txt"
rm $infoInput".tmpcol1.txt"
rm $infoInput".tmpcol2.txt"

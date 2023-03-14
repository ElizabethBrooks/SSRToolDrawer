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
cat $infoInput | tail -n+2 | sed '/^,/d' | sed 's/>//g' | sed 's/,/\t/g' | cut -f1-3,9 | sed 's/$/NEWLINE/g' > $infoOutput".tmp.txt"

# add the missing marker sequence info
cat $referenceInput | grep -A 1 "054-CM_015" | tail -1 >> $infoOutput".tmp.txt"

# combine the final two lines with the marker info
cat $infoOutput".tmp.txt" | sed 's/\n//g' | sed 's/NEWLINE/\n/g' > $infoOutput

# clean up
rm $infoOutput".tmp.txt"

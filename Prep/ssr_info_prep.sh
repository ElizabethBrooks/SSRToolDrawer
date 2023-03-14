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

# retrieve the missing sequence for marker 054-CM_015
markerSeqInfo=$(cat $referenceInput | grep -A 1 "054-CM_015" | tail -1)

# retrieve line with SSR info
markerSSRInfo=$(cat $infoInput | grep "054-CM_015")

# add the missing marker sequence info
# and retrieve the following fields for the SSR info csv file:
# Marker name, start repeat, end repeat, Sequence showing primer sequences and repeat in BOLD
# and exclude the header 
# and exclude the final three lines with notes by reversing the order of the file
# and fix 054-CM_015 marker tag by removing excess > symbol
# and convert the delimeter from commas to tabs
cat $infoInput | sed "s/$markerSSRInfo/$markerSSRInfo,$markerSeqInfo/g" | cut -d"," -f1-3,9 | tail -n+2 | tail -r | tail -n+4 | tail -r | sed 's/>//g' | sed 's/,/\t/g' > $infoOutput

#!/bin/bash

# script to convert SSR info csv files to tab delimeted text files
# usage: bash ssr_info_prep.sh inputsFile
# usage Ex: bash ssr_info_prep.sh inputPaths_ssr.txt

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve adapter absolute path for alignment
infoInput=$(grep "info:" ../"InputData/"$inputsFile | tr -d " " | sed "s/info://g")

# remove txt extension from input info file
infoOutput=$(echo $infoInput | sed 's/\.csv/\.txt/g')

# retrieve the following fields for the SSR info csv file:
# Marker name, start repeat, end repeat, Sequence showing primer sequences and repeat in BOLD
cat $infoInput | tail -n+2 | cut -d"," -f1-3,9 | sed 's/,/\t/g' > $infoOutput

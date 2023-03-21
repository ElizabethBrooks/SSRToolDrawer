#!/bin/bash

# script to convert SSR info csv files to tab delimeted text files
# usage: bash ssr_regions_prep.sh
# usage Ex: bash ssr_regions_prep.sh

# retrieve ssr info path
infoOutput=$(grep "info:" ../"InputData/inputPaths_ssr_pipeline.txt" | tr -d " " | sed "s/info://g")
#infoInput="/Users/bamflappy/GBCF/JRS/romero_Mar2023/Data/SSRinfo-2023.txt"
# retrieve ssr regions path
regionsPath=$(grep "ssrRegions:" ../"InputData/inputPaths_ssr_pipeline.txt" | tr -d " " | sed "s/ssrRegions://g")
#regionsPath="/Users/bamflappy/GBCF/JRS/romero_Mar2023/Data/SSRinfo-2023.bed"

# retrieve each marker, start, and end
cat $infoInput |  tr -d $'\r' | cut -f 1-3 > $regionsPath

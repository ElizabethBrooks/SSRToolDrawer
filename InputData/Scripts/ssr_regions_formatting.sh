#!/bin/bash

# script to convert SSR info csv files to tab delimeted text files
# usage: bash ssr_regions_formatting.sh
# usage Ex: bash ssr_regions_formatting.sh

# retrieve ssr info path
infoInput=$(grep "ssrInfo:" ../"InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrInfo://g")
#infoInput="/Users/bamflappy/GBCF/JRS/romero_Mar2023/Data/SSRinfo-2023.txt"
# retrieve ssr regions path
regionsPath=$(grep "ssrRegions:" ../"InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrRegions://g")
#regionsPath="/Users/bamflappy/GBCF/JRS/romero_Mar2023/Data/SSRinfo-2023.bed"

# remove txt extension from input primers file
inPath=$(echo $regionsPath | sed 's/\.bed//g')

# retrieve each marker, start, and end
cat $infoInput |  tr -d $'\r' | cut -f 1-3 > $inPath".tmp.txt"

# pre-clean up
rm $regionsPath

# loop over each line
while read line; do
	# retrieve marker name
	markerName=$(echo "$line" | cut -f 1)
	# retrieve start and end positions
	startPos=$(echo "$line" | cut -f 2)
	endPos=$(echo "$line" | cut -f 3)
	# subtract 1 from position for 0 indexed bed file
	startPos=$(($startPos - 1))
	endPos=$(($endPos -1 ))
	# output updated positions
	echo -e $markerName"\t"$startPos"\t"$endPos >> $regionsPath
done < $inPath".tmp.txt"

# clean up
rm $inPath".tmp.txt"

#!/bin/bash

# script to convert primer info csv files to tab delimeted text files
# usage: bash primer_info_formatting.sh
# usage Ex: bash primer_info_formatting.sh

# retrieve ssr primers path
outFile=$(grep "primers:" ../"inputs_ssr_pipeline.txt" | tr -d " " | sed "s/primers://g")
#outFile="/Users/bamflappy/GBCF/JRS/romero_Mar2023/Data/SSRprimer-info2-2023.bedpe"
# retrieve reference sequences path
referenceInput=$(grep "reference:" ../"inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")
#referenceInput="/Users/bamflappy/GBCF/JRS/romero_Mar2023/Data/finalset-2.fa"

# remove txt extension from input primers file
inPath=$(echo $outFile | sed 's/\.bedpe//g')

# set input primer sequence file
inFile=$inPath".csv"

# retrieve the missing marker sequence
markerSeq=$(cat $referenceInput | grep -A 1 "054-CM_015" | tail -1)

# retrieve the missing marker info
markerInfo=$(cat $inFile | grep "054-CM_015" | cut -d "," -f 1,5,6)

# retrieve the set of marker info fields
cat $inFile | tail -n+2 | sed '/^,/d' | grep -v "054-CM_015" | cut -d "," -f 1,5,6 > $inPath".tmpcol1.csv"

# retrieve the set of marker sequences
cat $inFile | tail -n+2 | sed '/^,/d' | grep -v "054-CM_015" | cut -d "," -f 8 > $inPath".tmpcol2.csv"

# add the missing marker info and sequence
paste -d "," $inPath".tmpcol1.csv" $inPath".tmpcol2.csv" > $inPath".tmp.csv"
echo  $markerInfo","$markerSeq >> $inPath".tmp.csv"

# pre-clean up
rm $outFile

# loop over each line
while read line; do
	# retrieve marker name
	markerName=$(echo $line | cut -d "," -f 1 | tr -d $'\r' | sed 's/^.//g')
	# retrieve the full sequence
	fullSeq=$(echo $line | cut -d "," -f 4 | tr -d $'\r')
	# retrieve the full sequence length
	fullLen=$(echo $fullSeq | awk '{ print length }')
	# retrieve sense primer sequence length
	senseLen=$(echo $line | cut -d "," -f 2 | tr -d $'\r' | awk '{ print length }')
	# retrieve antisense primer sequence length
	antiLen=$(echo $line | cut -d "," -f 3 | tr -d $'\r' | awk '{ print length }')
	# sense positions
	senseStart=0
	senseEnd=$(($senseLen-1))
	# antisense positions
	antiStart=$(($fullLen-$antiLen))
	antiStart=$(($antiStart-1))
	antiEnd=$(($fullLen-1))
	# output the sequence info in BEDPE format
	echo -e $markerName"\t"$senseStart"\t"$senseEnd"\t"$markerName"\t"$antiStart"\t"$antiEnd >> $outFile
done < $inPath".tmp.csv"

# clean up
rm $inPath".tmp.csv"
rm $inPath".tmpcol1.csv"
rm $inPath".tmpcol2.csv"

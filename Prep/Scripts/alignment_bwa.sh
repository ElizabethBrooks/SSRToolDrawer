#!/bin/bash

# script to perform bwa alignment of trimmed paired end reads
# usage: bash bwa_ssr_projects.sh inputsPath inputsFile baseDir

# required modules for ND CRC servers
#module load bio

# retrieve input argument of a inputs file
inputsPath=$1

# retrieve analysis outputs absolute path
inputsFile=$2

# retrieve base of working directory
baseDir=$3

# retrieve genome reference absolute path for alignment
ref=$(grep "reference:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")

# retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")

# retrieve analysis outputs path
outputsPath=$(dirname $inputsPath)

# name output file of inputs
versionFile=$outputsPath"/info/software_summary_prep.txt"

# output software version
echo "Alignment: " >> $versionFile
bwa &> tmp.txt
cat tmp.txt | head -3 >> $versionFile
rm tmp.txt
echo -e "\n" >> $versionFile

# make an outputs directory for analysis
anOut=$inputsPath"/aligned"

# this is used to re-start analyses
# check if the outputs directory already exsists
if [ -d $anOut ]; then
	echo "$anOut directory exists... exiting"
	exit
fi

# make the new directory
mkdir $anOut

# move to the outputs directory
#cd $anOut

# set trimmed reads absolute path
trimmedFolder=$inputsPath"/trimmed"

# loop through all forward and reverse paired reads and run bwa on each pair
# using 8 threads
for f1 in "$trimmedFolder"/*pForward.fq.gz; do
	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/.pForward\.fq\.gz//')
	# trim file path from current file name
	curSampleNoPath=$(basename $f1)
	curSampleNoPath=$(echo $curSampleNoPath | sed 's/.pForward\.fq\.gz//')
	# print status message
	echo "Processing $curSampleNoPath"
	# run bwa with default settings
	bwa mem -t 8 $ref $f1 $curSample"_pReverse.fq.gz" > $anOut"/"$curSampleNoPath".sam"
	# status message
	echo "Processed!"
done

# clean up
#rm -r "$trimmedFolder"

# print status message
echo "Analysis complete!"

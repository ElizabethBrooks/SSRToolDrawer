#!/bin/bash

# script to perform bwa alignment of trimmed paired end reads
# usage: bash bwa_ssr_projects.sh inputsFile inputsPath
# usage Ex: bash bwa_ssr_projects.sh inputPaths_romero_test_run1.txt

# required modules for ND CRC servers
#module load bio

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve input outputs path
inputsPath=$2

# retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
# retrieve genome reference absolute path for alignment
ref=$(grep "genomeReference:" ../"InputData/"$inputsFile | tr -d " " | sed "s/genomeReference://g")
# retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")

# name of output file of inputs
versionFile=$inputsPath"/software_prep_summary.txt"

# add software versions to outputs
bwa &> tmp.txt
cat tmp.txt | head -3 >> $versionFile
rm tmp.txt

# make an outputs directory for analysis
anOut=$inputsPath"/aligned"
mkdir $anOut

# move to the outputs directory
cd $anOut

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
	# create directory for current sample outputs
	mkdir "$curSampleNoPath"
	# print status message
	echo "Processing $curSampleNoPath"
	# run bwa with default settings
	bwa mem -t 8 $ref $f1 $curSample"_pReverse.fq.gz" > $curSampleNoPath".sam"
	# status message
	echo "Processed!"
done

# clean up
#rm -r "$trimmedFolder"

# print status message
echo "Analysis complete!"

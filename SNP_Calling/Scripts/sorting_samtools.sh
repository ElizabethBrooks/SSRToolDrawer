#!/bin/bash

# script to run the SSR pipeline
# usage: bash sorting_samtools.sh inputsPath baseDir

# retrieve inputs path
inputsPath=$1

# retrieve base of working directory
baseDir=$2

# retrieve analysis outputs path
outputsPath=$(dirname $inputsPath)

# name output file of inputs
versionFile=$outputsPath"/info/software_summary_SNP.txt"

# output software version
echo "Sorting: " >> $versionFile
samtools --version >> $versionFile
echo -e "\n" >> $versionFile

# set outputs path
outputsPath=$inputsPath"/sorted"
# create the directory
mkdir $outputsPath

# set inputs inputsPath
inputsPath=$inputsPath"/alignedFiltered50"


# Sorting Stage - SNP Calling Workflow

# loop through all filtered bam files
for f in $inputsPath"/"*".header.sam"; do
	# trim file path from current folder name
	curSampleNoPath=$(basename "$f" | sed 's/\.header\.sam$//g')
	# status message
	echo "Sorting $f"
	# convert sam to bam
	samtools view -@ 8 -bo $outputsPath"/"$curSampleNoPath".header.bam" $f
	# run samtools to prepare mapped reads for sorting
	samtools sort -@ 8 -n -o $outputsPath"/"$curSampleNoPath".sortedName.bam" -T "/tmp/"$curSampleNoPath".sortedName.bam" $outputsPath"/"$curSampleNoPath".header.bam"
	rm $outputsPath"/"$curSampleNoPath".header.bam"
	# run fixmate -m to update paired-end flags for singletons
	samtools fixmate -m $outputsPath"/"$curSampleNoPath".sortedName.bam" $outputsPath"/"$curSampleNoPath".sortedFixed.bam"
	rm $outputsPath"/"$curSampleNoPath".sortedName.bam"
	# run samtools to prepare mapped reads for sorting by coordinate
	samtools sort -@ 8 -o $outputsPath"/"$curSampleNoPath".sortedCoordinate.bam" -T "/tmp/"$curSampleNoPath".sortedCoordinate.bam" $outputsPath"/"$curSampleNoPath".sortedFixed.bam"
	rm $outputsPath"/"$curSampleNoPath".sortedFixed.bam"
	# remove duplicate reads
	#samtools markdup -r $outputsPath"/"$curSampleNoPath".sortedCoordinate.bam" $outputsPath"/"$curSampleNoPath".noDups.bam"
	#rm $outputsPath"/"$curSampleNoPath".sortedCoordinate.bam"
	# status message
	echo "Processed!"
done

# clean up
#rm -r $inputsPath

# status message
echo "Analysis conplete!"

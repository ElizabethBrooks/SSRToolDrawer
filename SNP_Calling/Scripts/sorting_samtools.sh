#!/bin/bash

# script to run the SSR pipeline
# usage: bash sorting_samtools.sh inputsPath projectDir

# retrieve inputs path
inputsPath=$1

# retrieve the project ID 
projectDir=$2

# name of output file of inputs
versionFile=$inputsPath"/software_VC_summary.txt"

# retrieve software version
samtools --version >> $versionFile


# Sorting Stage - SNP Calling Workflow

# set outputs path
outputsPath=$inputsPath"/sorted"
# create the directory
mkdir $outputsPath

# loop through all filtered bam files
for f in $inputsPath"/clipped/"*".readGroups.bam"; do
	# trim file path from current folder name
	curSampleNoPath=$(basename "$f" | sed 's/\.readGroups\.bam$//g')
	# status message
	echo "Sorting $f"
	# run samtools to prepare mapped reads for sorting
	samtools sort -@ 4 -n -o $outputsPath"/"$curSampleNoPath".sortedName.bam" -T "/tmp/"$curSampleNoPath".sortedName.bam" $f
	# run fixmate -m to update paired-end flags for singletons
	samtools fixmate -m $outputsPath"/"$curSampleNoPath".sortedName.bam" $outputsPath"/"$curSampleNoPath".sortedFixed.bam"
	#rm $outputsPath"/"$curSampleNoPath".sortedName.bam"
	# run samtools to prepare mapped reads for sorting by coordinate
	samtools sort -@ 4 -o $outputsPath"/"$curSampleNoPath".sortedCoordinate.bam" -T "/tmp/"$curSampleNoPath".sortedCoordinate.bam" $outputsPath"/"$curSampleNoPath".sortedFixed.bam"
	#rm $outputsPath"/"$curSampleNoPath".sortedFixed.bam"
	# remove duplicate reads
	samtools markdup -r $outputsPath"/"$curSampleNoPath".sortedCoordinate.bam" $outputsPath"/"$curSampleNoPath".noDups.bam"
	#rm $outputsPath"/"$curSampleNoPath".sortedCoordinate.bam"
	# status message
	echo "Processed!"
done

# status message
echo "Analysis conplete!"

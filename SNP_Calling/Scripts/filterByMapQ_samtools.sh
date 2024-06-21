#!/bin/bash

# script to run the SSR pipeline
# usage: bash filterByMapQ_samtools.sh inputsPath baseDir

# retrieve inputs path
inputsPath=$1

# retrieve base of working directory
baseDir=$2

# retrieve analysis outputs path
outputsPath=$(dirname $inputsPath)

# name output file of inputs
versionFile=$outputsPath"/info/software_summary_SNP.txt"

# output software version
echo "Alignment filtering: " >> $versionFile
samtools --version >> $versionFile

# set outputs path
outputsPath=$inputsPath"/alignedFilteredMapQ"
# create the directory
mkdir $outputsPath

# set inputs inputsPath
inputsPath=$inputsPath"/sorted"


# Alignments Filtering Stage - SNP Calling Workflow

# status message
echo "Performing alignment filtering..."

# loop over each coordinate sorted sample
for f in $inputsPath"/"*".sortedCoordinate.bam"; do
	# trim file path from current folder name
	curSampleNoPath=$(basename "$f" | sed 's/\.sortedCoordinate\.bam$//g')
	# status message
	echo "Filtering $f"
	# keep only unique read alignments using a mapq score of 60
	samtools view -@ 8 -bq 60 $f > $outputsPath"/"$curSampleNoPath".filteredMapQ.bam"
	# index bam file
	samtools index -@ 8 $outputsPath"/"$curSampleNoPath".filteredMapQ.bam"
done

# clean up
#rm -r $inputsPath

# status message
echo "Analysis conplete!"

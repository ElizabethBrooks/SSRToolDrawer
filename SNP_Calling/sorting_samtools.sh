#!/bin/bash

# script to run the SSR pipeline
# usage: bash sorting_samtools.sh inputsFile inputsPath

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve input outputs path
inputsPath=$2

# retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")

# setup the variant calling directory
dataPath=$inputsPath"/sorted"
# create the directory
mkdir $dataPath

# name of output file of inputs
versionFile=$inputsPath"/software_VC_summary.txt"
# add pipeline info to outputs
echo -e "SSR pipeline variant calling software versions for $projectDir \n" > $versionFile

# retrieve software version
samtools --version >> $versionFile


# Sorting Stage - SNP Calling Workflow

# set input paths
inputsPath=$inputsPath"/aligned"

# loop through all filtered sam files
for f in $inputsPath"/"*".header.filter50.sam"; do
	# remove the file extension
	path=$(echo $f | sed 's/\.sam$//g')
	# trim file path from current folder name
	curSampleNoPath=$(basename "$f" | sed 's/\.sam$//g')
	# status message
	echo "Processing file "$path".sam ..."
	# convert output sam files to bam format for downstream analysis
	samtools view -@ 4 -bS $f > $path".bam"
	# run samtools to prepare mapped reads for sorting
	# using 8 threads
	samtools sort -@ 4 -n -o $path".sortedName.bam" -T "/tmp/"$curSampleNoPath".sortedName.bam" $path".bam"
	# run fixmate -m to update paired-end flags for singletons
	samtools fixmate -m $path".sortedName.bam" $path".sortedFixed.bam"
	# run samtools to prepare mapped reads for sorting by coordinate
	# using 8 threads
	samtools sort -@ 4 -o $path".sortedCoordinate.bam" -T "/tmp/"$curSampleNoPath".sortedCoordinate.bam" $path".sortedFixed.bam"
	# remove duplicate reads
	samtools markdup -r $path".sortedCoordinate.bam" $path".noDups.bam"
	# clean up
	rm $f
	rm $path".sortedName.bam"
	rm $path".sortedFixed.bam"
	rm $path".sortedCoordinate.bam"
	# status message
	echo "Processed!"
done

# status message
echo "Analysis conplete!"

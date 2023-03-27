#!/bin/bash

# script to clip primer and ssr sequences
# usage: bash clipping_samtools_bamclipper.sh inputsPath baseDir

# retrieve inputs path
inputsPath=$1

# retrieve base of working directory
baseDir=$2

# retrieve primers path
primerPath=$(grep "primers:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/primers://g")
# bamclipper tool path
clipperPath=$(grep "bamclipperTool:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/bamclipperTool://g")

# name of output file of inputs
versionFile=$inputsPath"/software_VC_summary.txt"

# retrieve software version
samtools --version >> $versionFile


# Clipping Stage - SNP Calling Workflow

# set outputs path
outputsPath=$inputsPath"/clipped"
# create the directory
mkdir $outputsPath

# copy bamclipper software directory
cp -r $clipperPath $outputsPath

# retrieve bamclipper directory name
clipperBase=$(basename $clipperPath)

# move to outputs directory
cd $outputsPath"/"$clipperBase

# loop through all aligned sam files
for f1 in $inputsPath"/sorted/"*".sortedCoordinate.bam"; do
	# trim file path from current folder name
	curSampleNoPath=$(basename "$f1" | sed 's/\.sortedCoordinate\.bam$//g')
	# print status message
	echo "Clipping $f1"
	# index the bam file
	samtools index $f1 
	# soft mask primers sequences
	./bamclipper.sh -b $f1 -p $primerPath -n 4
	# add read groups
	samtools addreplacerg -@ 4 -r ID:"SSR_"$runNum"_"$curSampleNoPath -r SM:$curSampleNoPath -o $outputsPath"/"$curSampleNoPath".readGroups.bam" $outputsPath"/"$clipperBase"/"$curSampleNoPath".sortedCoordinate.primerclipped.bam"
	# status message
	echo "Processed!"
done

# clean up
rm -r $outputsPath"/"$clipperBase

# status message
echo "Analysis complete!"

#!/bin/bash

# script to clip primer and ssr sequences
# usage: bash clipping_samtools_bamclipper.sh inputsPath

# retrieve inputs path
inputsPath=$1

# name of output file of inputs
versionFile=$inputsPath"/software_VC_summary.txt"

# retrieve software version
samtools --version >> $versionFile


# Clipping Stage - SNP Calling Workflow

# set outputs path
outputsPath=$inputsPath"/clipped"

# loop through all aligned sam files
for f1 in $inputsPath"/filtered/"*".header.sam"; do
	# trim file path from current folder name
	curSampleNoPath=$(basename "$f1" | sed 's/\.header\.sam$//g')
	# print status message
	echo "Processing $f1"
	# convert sam to bam
	#samtools view -@ 4 -bo $outputsPath"/"$curSampleNoPath".header.bam" $f1
	# coordinate sort the filtered sequences
	#samtools sort -@ 4 -o $outputsPath"/"$curSampleNoPath".sortedCoordinate.bam" -T "/tmp/"$curSampleNoPath".sortedCoordinate.bam" $outputsPath"/"$curSampleNoPath".header.bam"
	#rm $outputsPath"/"$curSampleNoPath".header.bam"
	# index the bam file
	#samtools index $outputsPath"/"$curSampleNoPath".sortedCoordinate.bam" 
	# remove primers sequences
	./bamclipper.sh -b $outputsPath"/"$curSampleNoPath".sortedCoordinate.bam" -p $primerPath -n 4
	#rm $outputsPath"/"$curSampleNoPath".sortedCoordinate.bam"
	# remove SSR sequences
	samtools view -@ 4 -bo $outputsPath"/"$curSampleNoPath".overlap.bam" -U $outputsPath"/"$curSampleNoPath".noOverlap.bam" -L $regionsPath $outputsPath"/"$curSampleNoPath".sortedCoordinate.primerclipped.bam"
	#rm $outputsPath"/"$curSampleNoPath".sortedCoordinate.primerclipped.bam"
	#rm $outputsPath"/"$curSampleNoPath".overlap.bam"
	# add read groups
	samtools addreplacerg -@ 4 -r ID:"SSR_"$runNum"_"$curSampleNoPath -r SM:$curSampleNoPath -o $outputsPath"/"$curSampleNoPath".readGroups.bam" $outputsPath"/"$curSampleNoPath".noOverlap.bam"
	#rm $outputsPath"/"$curSampleNoPath".noOverlap.bam"
	# status message
	echo "Processed!"
done

# status message
echo "Analysis complete!"

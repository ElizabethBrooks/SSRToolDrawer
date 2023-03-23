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

# copy pipeline scripts to the clipped directory
cp -r $clipperPath"/"* $outputsPath

# move to outputs path
cd $outputsPath

# loop through all aligned sam files
for f1 in $inputsPath"/filtered/"*".noDups.bam"; do
	# trim file path from current folder name
	curSampleNoPath=$(basename "$f1" | sed 's/\.noDups\.bam$//g')
	# print status message
	echo "Processing $f1"
	# index the bam file
	samtools index $outputsPath"/"$curSampleNoPath".noDups.bam" 
	# soft mask primers sequences
	./bamclipper.sh -b $outputsPath"/"$curSampleNoPath".noDups.bam" -p $primerPath -n 4
	#rm $outputsPath"/"$curSampleNoPath".noDups.bam"
	# remove SSR sequences
	#samtools view -@ 4 -bo $outputsPath"/"$curSampleNoPath".overlap.bam" -U $outputsPath"/"$curSampleNoPath".noOverlap.bam" -L $regionsPath $outputsPath"/"$curSampleNoPath".noDups.primerclipped.bam"
	#rm $outputsPath"/"$curSampleNoPath".noDups.primerclipped.bam"
	#rm $outputsPath"/"$curSampleNoPath".overlap.bam"
	# add read groups
	samtools addreplacerg -@ 4 -r ID:"SSR_"$runNum"_"$curSampleNoPath -r SM:$curSampleNoPath -o $outputsPath"/"$curSampleNoPath".readGroups.bam" $outputsPath"/"$curSampleNoPath".noDups.primerclipped.bam"
	# status message
	echo "Processed!"
done

# status message
echo "Analysis complete!"

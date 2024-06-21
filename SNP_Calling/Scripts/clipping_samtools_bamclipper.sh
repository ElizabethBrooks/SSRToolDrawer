#!/bin/bash

# script to clip primer and ssr sequences
# usage: bash clipping_samtools_bamclipper.sh inputsPath baseDir runNum

# retrieve inputs path
inputsPath=$1

# retrieve base of working directory
baseDir=$2

# retrieve run number
sampleRun=$3

# retrieve primers path
primerPath=$(grep "primers:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/primers://g")
# bamclipper tool path
clipperPath=$(grep "bamclipperTool:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/bamclipperTool://g")
# retrieve analysis outputs path
infoPath=$(grep "outputs:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/outputs://g")

# name output file of inputs
versionFile=$infoPath"/info/software_summary_SNP.txt"

# output software version
echo "Clipping: " >> $versionFile
samtools --version >> $versionFile
echo -e "\n" >> $versionFile

# set outputs path
outputsPath=$inputsPath"/clipped"
# create the directory
mkdir $outputsPath

# set inputs inputsPath
inputsPath=$inputsPath"/alignedFilteredMapQ"


# Clipping Stage - SNP Calling Workflow

# add run tags to sample files names
for i in $inputsPath"/"*".filteredMapQ.bam"; do
	# retrieve sanmepl run tag
	newName=$(echo $i | sed "s/\.filteredMapQ\.bam$/\_$sampleRun\.filteredMapQ\.bam/g")
	# update the sample file name
	mv $i $newName
done

# copy bamclipper software directory
cp -r $clipperPath $outputsPath

# retrieve bamclipper directory name
clipperBase=$(basename $clipperPath)

# move to outputs directory
cd $outputsPath"/"$clipperBase

# loop through all aligned sam files
for f1 in $inputsPath"/"*".filteredMapQ.bam"; do
	# trim file path from current folder name
	curSampleNoPath=$(basename "$f1" | sed 's/\.filteredMapQ\.bam$//g')
	# print status message
	echo "Clipping $f1"
	# index the bam file
	samtools index $f1 
	# soft mask primers sequences
	./bamclipper.sh -b $f1 -p $primerPath -n 4
	# add read groups
	samtools addreplacerg -@ 8 -r ID:"SSR_"$runNum"_"$curSampleNoPath -r SM:$curSampleNoPath -o $outputsPath"/"$curSampleNoPath".readGroups.bam" $outputsPath"/"$clipperBase"/"$curSampleNoPath".filteredMapQ.primerclipped.bam"
	# status message
	echo "Processed!"
done

# clean up
#rm -r $inputsPath
rm -r $outputsPath"/"$clipperBase

# move back to pipeline scripts directory
cd $baseDir"/SNP_Calling/Scripts"

# status message
echo "Analysis complete!"

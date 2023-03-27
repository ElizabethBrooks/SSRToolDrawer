#!/bin/bash

# script to perform trimming of paired end reads
# usage: bash trimmomatic_ssr_projects.sh inputsFile outputsPath

# required modules for ND CRC servers
#module load bio

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve input outputs path
outputsPath=$2

# retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")
# retrieve adapter absolute path for alignment
adapterPath=$(grep "adapter:" ../../"InputData/inputPaths_ssr_pipeline.txt" | tr -d " " | sed "s/adapter://g")

# name of output file of inputs
versionFile=$outputsPath"/software_prep_summary.txt"

# add software version to outputs
echo "Trimmomatic:" >> $versionFile
trimmomatic -version >> $versionFile

# make a new directory for analysis
trimOut=$outputsPath"/trimmed"
mkdir $trimOut

# move to the new directory
#cd $trimOut

# set directory for qc input
qcFolder=$outputsPath"/qc"

# loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*"_R1_001.fastq.gz"; do
	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R._001\.fastq\.gz//')
	# set paired file name
	f2=$curSample"_R2_001.fastq.gz"
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_R._001\.fastq\.gz//')
	# print status message
	echo "Processing $sampleTag"
	# determine phred score for trimming
	if grep -iF "Illumina 1.5" $qcFolder"/"$sampleTag"_R1_001_fastqc/fastqc_data.txt"; then
		score=64
	elif grep -iF "Illumina 1.9" $qcFolder"/"$sampleTag"_R1_001_fastqc/fastqc_data.txt"; then
		score=33
	else
		echo "ERROR: Illumina encoding not found... exiting"
		exit 1
	fi
	echo $qcFolder"/"$sampleTag"_R1_001_fastqc/fastqc_data.txt"
	# perform adapter trimming on paired reads using 8 threads
	trimmomatic PE -threads 4 -phred"$score" $f1 $f2 $trimOut"/"$sampleTag"_pForward.fq.gz" $trimOut"/"$sampleTag"_uForward.fq.gz" $trimOut"/"$sampleTag"_pReverse.fq.gz" $trimOut"/"$sampleTag"_uReverse.fq.gz" ILLUMINACLIP:"$adapterPath" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	# print status message
	echo "Processed!"
done

# clean up
#rm -r "$trimmedFolder"

# print status message
echo "Analysis complete!"

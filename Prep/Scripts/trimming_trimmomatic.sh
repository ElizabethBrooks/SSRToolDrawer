#!/bin/bash

# script to perform trimming of paired end reads
# usage: bash trimmomatic_ssr_projects.sh inputsPath inputsFile baseDir

# required modules for ND CRC servers
#module load bio

# retrieve input argument of a inputs file
inputsPath=$1

# retrieve analysis outputs absolute path
inputsFile=$2

# retrieve base of working directory
baseDir=$3

# retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")

# retrieve adapter absolute path for alignment
adapterPath=$(grep "adapter:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/adapter://g")

# retrieve analysis outputs path
outputsPath=$(dirname $inputsPath)

# name output file of inputs
versionFile=$outputsPath"/info/software_summary_prep.txt"

# output software version
echo "Read trimming: " >> $versionFile
trimmomatic -version >> $versionFile
echo -e "\n" >> $versionFile

# make a new directory for analysis
trimOut=$inputsPath"/trimmed"
mkdir $trimOut

# move to the new directory
#cd $trimOut

# set directory for qc input
qcFolder=$inputsPath"/qc"

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
	trimmomatic PE -threads 8 -phred"$score" $f1 $f2 $trimOut"/"$sampleTag"_pForward.fq.gz" $trimOut"/"$sampleTag"_uForward.fq.gz" $trimOut"/"$sampleTag"_pReverse.fq.gz" $trimOut"/"$sampleTag"_uReverse.fq.gz" ILLUMINACLIP:"$adapterPath" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	# print status message
	echo "Processed!"
done

# clean up
#rm -r "$trimmedFolder"

# print status message
echo "Analysis complete!"

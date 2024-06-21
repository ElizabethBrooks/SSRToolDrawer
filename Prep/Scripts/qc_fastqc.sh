#!/bin/bash

# script to perform fastqc quality control of paired end reads
# usage: bash fastqc_ssr_projects.sh inputsPath inputsFile baseDir

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

# retrieve analysis outputs path
outputsPath=$(dirname $inputsPath)

# name output file of inputs
versionFile=$outputsPath"/info/software_summary_prep.txt"

# output software version
echo "QC: " >> $versionFile
fastqc -version >> $versionFile
echo -e "\n" >> $versionFile

# make a new directory for analysis
qcOut=$inputsPath"/qc"
mkdir $qcOut

# move to the new directory
#cd $qcOut

# loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*"_R1_001.fastq.gz"; do
	# trim path from file name
	noPath=$(basename $f1 | sed 's/_R._001\.fastq\.gz//')
	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R._001\.fastq\.gz//')
	# set paired file name
	f2=$curSample"_R2_001.fastq.gz"
	# print status message
	echo "Processing $noPath"
	# perform QC on both paired end reads for the current sample
	fastqc $f1 -o $qcOut --extract
	fastqc $f2 -o $qcOut --extract
	# print status message
	echo "Processed!"
done

# print status message
echo "Analysis complete!"

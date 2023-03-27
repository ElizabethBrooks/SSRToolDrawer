#!/bin/bash

# script to perform fastqc quality control of paired end reads
# usage: bash fastqc_ssr_projects.sh inputsFile outputsPath

# required modules for ND CRC servers
# module load bio

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve input outputs path
outputsPath=$2

# retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")

# name of output file of inputs
versionFile=$outputsPath"/software_prep_summary.txt"

# report software version
fastqc -version >> $versionFile

# make a new directory for analysis
qcOut=$outputsPath"/qc"
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

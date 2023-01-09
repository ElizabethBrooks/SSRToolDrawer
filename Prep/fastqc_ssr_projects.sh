#!/bin/bash
#Script to perform fastqc quality control of paired end reads
#Usage: bash fastqc_ssr_projects.sh inputsFile
#Usage Ex: bash fastqc_ssr_projects.sh inputPaths_romero_run1.txt

#Required modules for ND CRC servers
#module load bio

#Retrieve input argument of a inputs file
inputsFile=$1

#Retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
#Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

#Directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_workflow_prep"

#Name of output file of inputs
inputOutFile=$outputsPath"/pipeline_prep_summary.txt"
versionFile=$outputsPath"/software_prep_summary.txt"

#Report software version
fastqc -version >> $versionFile

#Make a new directory for analysis
qcOut=$outputsPath"/qc"
mkdir $qcOut

#Move to the new directory
cd $qcOut

#Loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*_R1_001.fastq.gz; do
	#Trim path from file name
	noPath=$(basename $f1 | sed 's/_R._001\.fastq\.gz//')
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R._001\.fastq\.gz//')
	#Set paired file name
	f2=$curSample"_R2_001.fastq.gz"
	#Print status message
	echo "Processing $noPath"
	#Perform QC on both paired end reads for the current sample
	fastqc $f1 -o $qcOut --extract
	fastqc $f2 -o $qcOut --extract
	#Output run inputs
	echo "fastqc $f1 -o $qcOut --extract" >> $inputOutFile
	echo "fastqc $f2 -o $qcOut --extract" >> $inputOutFile
	#Print status message
	echo "Processed!"
done

#Print status message
echo "Analysis complete!"

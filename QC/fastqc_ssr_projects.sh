#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N fastqc_ssr_projects_jobOutput
#$ -pe smp 8
#Script to perform fastqc quality control of paired end reads
#Usage: qsub fastqc_ssr_projects.sh inputsFile
#Usage Ex: qsub fastqc_ssr_projects.sh inputPaths_romero_test_July2022.txt

#Required modules for ND CRC servers
module load bio

#Retrieve input argument of a inputs file
inputsFile=$1

#Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")
#Retrieve adapter absolute path for alignment
adapterPath=$(grep "adapter:" ../"InputData/"$inputsFile | tr -d " " | sed "s/adapter://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

#Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir
mkdir $outputsPath

#Make a new directory for analysis
qcOut=$outputsPath"/qc"
mkdir $qcOut
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $qcOut directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to the new directory
cd $qcOut

#Report software version
fastqc -version > "summary.txt"
#Loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*_R1_001.fastq.gz; do
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R._001\.fastq\.gz//')
	#Set paired file name
	f2=$curSample"_R2_001.fastq.gz"
	#Print status message
	echo "Processing $curSample"
	#Perform QC on both paired end reads for the current sample
	fastqc $f1 -o $qcOut --extract
	fastqc $f2 -o $qcOut --extract
	#Print status message
	echo "Processed!"
done
#Print status message
echo "Analysis complete!"

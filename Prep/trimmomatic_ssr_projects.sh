#!/bin/bash
#Script to perform trimming of paired end reads
#Usage: bash trimmomatic_ssr_projects.sh inputsFile
#Usage Ex: bash trimmomatic_ssr_projects.sh inputPaths_romero_test_run1.txt

#Required modules for ND CRC servers
#module load bio

#Retrieve input argument of a inputs file
inputsFile=$1

#Retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
#Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")
#Retrieve adapter absolute path for alignment
adapterPath=$(grep "adapter:" ../"InputData/"$inputsFile | tr -d " " | sed "s/adapter://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

#Directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_basicWorkflow"

#Name of output file of inputs
inputOutFile=$outputsPath"/pipeline_summary.txt"
versionFile=$outputsPath"/software_summary.txt"
#Add software version to outputs
trimmomatic -version >> $versionFile

#Make a new directory for analysis
trimOut=$outputsPath"/trimmed"
mkdir $trimOut

#Move to the new directory
cd $trimOut

#Loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*_R1_001.fastq.gz; do
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R._001\.fastq\.gz//')
	#Set paired file name
	f2=$curSample"_R2_001.fastq.gz"
	#Trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_R._001\.fastq\.gz//')
	#Print status message
	echo "Processing $sampleTag"
	#Determine phred score for trimming
	if grep -iF "Illumina 1.5" $outputsPath"/qc/"$sampleTag"_R1_001_fastqc/fastqc_data.txt"; then
		score=64
	elif grep -iF "Illumina 1.9" $outputsPath"/qc/"$sampleTag"_R1_001_fastqc/fastqc_data.txt"; then
		score=33
	else
		echo "ERROR: Illumina encoding not found... exiting"
		#echo "ERROR: Illumina encoding not found for $curSample" >> $inputOutFile
		exit 1
	fi
	echo $outputsPath"/qc/"$sampleTag"_R1_001_fastqc/fastqc_data.txt"
	#Perform adapter trimming on paired reads using 8 threads
	trimmomatic PE -threads 8 -phred"$score" $f1 $f2 $sampleTag"_pForward.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag"_pReverse.fq.gz" $sampleTag"_uReverse.fq.gz" ILLUMINACLIP:"$adapterPath" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	#Add run inputs to output summary file
	echo trimmomatic PE -threads 8 -phred"$score" $f1 $f2 $sampleTag"_pForward.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag"_pReverse.fq.gz" $sampleTag"_uReverse.fq.gz" ILLUMINACLIP:"$adapterPath" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >> $inputOutFile
	#Clean up
	rm -r $noPath"_R1_001_fastqc.zip"
	rm -r $noPath"_R1_001_fastqc/"
	rm -r $noPath"_R2_001_fastqc.zip"
	rm -r $noPath"_R2_001_fastqc/"
	#Print status message
	echo "Processed!"
done

#Print status message
echo "Analysis complete!"

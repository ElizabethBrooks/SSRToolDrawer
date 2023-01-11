#!/bin/bash
#Script to perform bwa alignment of trimmed paired end reads
#Usage: bash bwa_ssr_projects.sh inputsFile
#Usage Ex: bash bwa_ssr_projects.sh inputPaths_romero_test_run1.txt

#Required modules for ND CRC servers
#module load bio

#Retrieve input argument of a inputs file
inputsFile=$1

# retrieve input outputs path
outputsPath=$2

#Retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
#Retrieve genome reference absolute path for alignment
ref=$(grep "genomeReference:" ../"InputData/"$inputsFile | tr -d " " | sed "s/genomeReference://g")
#Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")

#Name of output file of inputs
versionFile=$outputsPath"/software_prep_summary.txt"

#Add software versions to outputs
bwa &> tmp.txt
cat tmp.txt | head -3 >> $versionFile
rm tmp.txt

#Make an outputs directory for analysis
anOut=$outputsPath"/aligned"
mkdir $anOut

#Move to the outputs directory
cd $anOut

#Set trimmed reads absolute path
trimmedFolder=$outputsPath"/trimmed"

#Loop through all forward and reverse paired reads and run bwa on each pair
# using 8 threads
for f1 in "$trimmedFolder"/*pForward.fq.gz; do
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/.pForward\.fq\.gz//')
	#Trim file path from current file name
	curSampleNoPath=$(basename $f1)
	curSampleNoPath=$(echo $curSampleNoPath | sed 's/.pForward\.fq\.gz//')
	#Create directory for current sample outputs
	mkdir "$curSampleNoPath"
	#Print status message
	echo "Processing $curSampleNoPath"
	#Run bwa with default settings
	bwa mem -t 8 $ref $f1 $curSample"_pReverse.fq.gz" > $curSampleNoPath".sam"
	# status message
	echo "Processed!"
done

#Clean up
rm -r "$trimmedFolder"

#Print status message
echo "Analysis complete!"

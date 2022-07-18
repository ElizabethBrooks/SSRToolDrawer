#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_pipeline_jobOutput
#$ -pe smp 8
#Script to run the SSR pipeline
#Usage: qsub ssr_pipeline.sh inputsFile
#Usage Ex: qsub ssr_pipeline.sh inputPaths_romero_July2022.txt

#Required modules for ND CRC servers
module load bio

#Activate the python2 environment for local run
source /afs/crc.nd.edu/user/e/ebrooks5/.bashrc
conda activate /afs/crc.nd.edu/user/e/ebrooks5/.conda/envs/python2

#Activate the python2 environment for job run
#conda activate python2

#Make sure numpy is installed
#pip install numpy

#Retrieve input argument of a inputs file
inputsFile=$1

#Retrieve paired reads absolute path for alignment
readPath=$(grep "reads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/reads://g")
#Retrieve adapter absolute path for alignment
infoPath=$(grep "info:" ../"InputData/"$inputsFile | tr -d " " | sed "s/info://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

#Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir

#Set input paths
inputsPath=$outputsPath"/aligned"

#Make an outputs directory for analysis
anOut=$outputsPath"/SSR"
mkdir $anOut
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $anOut directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to the outputs directory
cd $anOut

#Name output file of inputs
inputOutFile="summary.txt"
#Add pipeline info to outputs
echo "SSR pipline for $projectDir" > $inputOutFile

#Loop through all filtered sam files
for f1 in "$inputsPath"/*_L001.sam; do
	#Print status message
	echo "Processing $f1"
	#Run SSR pipeline
	python2 GapGenes.v3.py -sam $f1 -C $infoPath -P "unpaired"
	python2 SnipMatrix.py $f1".Matrix.txt"
	#Write inputs out to summary file
	echo python2 GapGenes.v3.py -sam $f1 -C $infoPath -P "unpaired" >> $inputOutFile
	echo python2 SnipMatrix.py $f1".Matrix.txt" >> $inputOutFile
	#Status message
	echo "Processed!"
done

#Retrieve sample list
#sampleList=$(ls /scratch365/ebrooks5/romero_test_July2022/sam/*_L001.sam.filter50.sam | sed "s/^/\"/g" | sed "s/\.sam\.filter50\.sam/\",/g" | tr '\n' ' ')
#sampleList=$(ls "$readPath"/*_L001.sam.filter50.sam | sed "s/\/scratch365\/ebrooks5\/romero_test_July2022\/sam\//\"/g" | sed "s/\.sam\.filter50\.sam/\",/g" | tr '\n' ' ')

#Format matrix
#python2 Format_Matrix.py
#echo python2 Format_Matrix.py >> $inputOutFile

#Print status message
#echo "Analysis complete!"

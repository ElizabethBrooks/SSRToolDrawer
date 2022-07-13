#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_pipeline_jobOutput
#$ -pe smp 8
#Script to perform trimmomatic trimming of paired end reads
#Usage: qsub ssr_pipeline.sh inputsFile
#Usage Ex: qsub ssr_pipeline.sh inputPaths_romero_test_July2022.txt

#Required modules for ND CRC servers
module load bio

#Activate the python2 environment
conda activate python2

#Retrieve input argument of a inputs file
inputsFile=$1

#Retrieve paired reads absolute path for alignment
readPath=$(grep "reads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/reads://g")
#Retrieve adapter absolute path for alignment
infoPath=$(grep "info:" ../"InputData/"$inputsFile | tr -d " " | sed "s/info://g")

#Make a new directory for project analysis
projectDir=$(basename $readPath)

#Move to the data directory
cd $readPath

#Name output file of inputs
inputOutFile="summary.txt"
echo "SSR pipline for $projectDir" > $inputOutFile

#Loop through all filtered sam files
for f1 in "$readPath"/*_L001.sam.filter50.sam; do
	#Print status message
	echo "Processing $curSample"
	#Run SSR pipeline
	python GapGenes.v3.py -sam $f1 -C $infoPath -P "unpaired"
	python SnipMatrix.py $f1".Matrix.txt"
	python Format_Matrix.py
	#Write inputs out to summary file
	echo python GapGenes.v3.py -sam $f1 -C $infoPath -P "unpaired" >> $inputOutFile
	echo python SnipMatrix.py $f1".Matrix.txt" >> $inputOutFile
	echo python Format_Matrix.py >> $inputOutFile
	#Status message
	echo "Processed!"
done
#Print status message
echo "Analysis complete!"

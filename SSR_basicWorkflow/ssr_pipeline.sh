#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_pipeline_jobOutput
#Script to run the SSR pipeline
#Usage: qsub ssr_pipeline.sh inputsFile
#Usage Ex: qsub ssr_pipeline.sh inputPaths_romero_run1.txt

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

#Prepare reads for SSR analysis
cd ../Prep
bash fastqc_ssr_projects.sh $1
bash trimmomatic_ssr_projects.sh $1
bash bwa_ssr_projects.sh $1

#Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")
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

#Copy pipeline scripts to outputs directory
cd ../SSR_basicWorkflow
cp GapGenes.v3.py $inputsPath
cp SnipMatrix.py $inputsPath
cp Format_Matrix.py $inputsPath

#Move to the outputs directory
cd $inputsPath

#Name output file of inputs
inputOutFile=$anOut"/summary.txt"
#Add pipeline info to outputs
echo "SSR pipline for $projectDir" > $inputOutFile

#Loop through all filtered sam files
for f1 in "$inputsPath"/*.sam; do
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
#for i in "$inputsPath"/*.sam; do basename $i | sed "s/^/\"/g" | sed "s/\.sam/\",/g" | tr '\n' ' '; done

#Format matrix
python2 Format_Matrix.py
echo python2 Format_Matrix.py >> $inputOutFile

#Clean up
rm "$inputsPath"/*.py

#Print status message
echo "Analysis complete!"

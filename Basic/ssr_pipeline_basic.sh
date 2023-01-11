#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_basic_jobOutput
#Script to run the SSR pipeline
#Usage: qsub ssr_pipeline_basic.sh inputsFile
#Usage Ex: qsub ssr_pipeline_basic.sh inputPaths_romero_run1.txt
#Usage Ex: qsub ssr_pipeline_basic.sh inputPaths_romero_run2.txt
#Usage Ex: qsub ssr_pipeline_basic.sh inputPaths_romero_run3.txt
#Usage Ex: qsub ssr_pipeline_basic.sh inputPaths_romero_run4.txt
#Usage Ex: qsub ssr_pipeline_basic.sh inputPaths_romero_run5.txt

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

#Retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
#Retrieve adapter absolute path for alignment
infoPath=$(grep "info:" ../"InputData/"$inputsFile | tr -d " " | sed "s/info://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

#Make a new directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_basic"
mkdir $outputsPath
#Check if the folder already exists
#if [ $? -ne 0 ]; then
#	echo "The $outputsPath directory already exsists... please remove before proceeding."
#	exit 1
#fi

# setup the inputs path
inputsPath=$outputsPath"/"$projectDir"_SSR_prep"
mkdir $inputsPath
#Check if the folder already exists
#if [ $? -ne 0 ]; then
#	echo "The $inputsPath directory already exsists... please remove before proceeding."
#	exit 1
#fi

# prepare data for analysis
#cd ../Prep
#bash ssr_pipeline_prep.sh $inputsFile $inputsPath


#SSR Analysis Stage - Basic Workflow

# status message
echo "SSR basic analysis started..."

#Set input paths
inputsPath=$inputsPath"/aligned"

#Copy pipeline scripts to inputs directory
cd ../Basic
cp GapGenes.v3.py $inputsPath
cp SnipMatrix.py $inputsPath
cp Format_Matrix.py $inputsPath

#Move to the inputs directory
cd $inputsPath

#Loop through all filtered sam files
for f1 in "$inputsPath"/*.sam; do
	#Print status message
	echo "Processing $f1"
	#Run SSR pipeline
	python2 GapGenes.v3.py -sam $f1 -C $infoPath -P "unpaired"
	python2 SnipMatrix.py $f1".Matrix.txt"
	#Status message
	echo "Processed!"
done

#Retrieve and format sample tag list
sampleTags=$(for i in "$inputsPath"/*.sam; do basename $i | sed "s/^/\"/g" | sed "s/\.sam/\",/g" | tr '\n' ' '; done)
sampleTags=$(echo $sampleTags | sed 's/.$//')

#Find and replace the sample list
sed -i "s/\"FIND_ME_REPLACE_ME\"/$sampleTags/g" Format_Matrix.py

#Format matrix
python2 Format_Matrix.py

#Re-name and move output matrix
mv SNP_Matrix.txt $outputsPath"/"$projectDir"_SNP_Matrix.txt"

# clean up
#inputsPath=$(dirname $inputsPath)
#rm -r $inputsPath

# status message
echo "SSR SNP analysis complete!"

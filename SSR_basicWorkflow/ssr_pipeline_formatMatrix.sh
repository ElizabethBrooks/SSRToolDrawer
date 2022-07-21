#!/bin/bash
#Script to run the SSR pipeline
#Usage: bash ssr_pipeline_formatMatrix.sh inputsFile
#Usage Ex: bash ssr_pipeline_formatMatrix.sh inputPaths_romero_run1.txt
#Usage Ex: bash ssr_pipeline_formatMatrix.sh inputPaths_romero_run2.txt
#Usage Ex: bash ssr_pipeline_formatMatrix.sh inputPaths_romero_run3.txt
#Usage Ex: bash ssr_pipeline_formatMatrix.sh inputPaths_romero_run4.txt
#Usage Ex: bash ssr_pipeline_formatMatrix.sh inputPaths_romero_run5.txt

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
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

#Make a new directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_basicWorkflow"
mkdir $outputsPath
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi


#SSR Analysis Stage

#Set input paths
inputsPath=$outputsPath"/aligned"

#Move to directory with SSR analysis scripts
cd ../SSR_basicWorkflow

#Copy pipeline scripts to inputs directory
cp Format_Matrix.py $inputsPath

#Move to the inputs directory
cd $inputsPath

#Retrieve and format sample tag list
sampleTags=$(for i in "$inputsPath"/*.sam; do basename $i | sed "s/^/\"/g" | sed "s/\.sam/\",/g" | tr '\n' ' '; done)
sampleTags=$(echo $sampleTags | sed 's/.$//')

#Find and replace the sample list
sed -i "s/\"FIND_ME_REPLACE_ME\"/$sampleTags/g" Format_Matrix.py

#Format matrix
python2 Format_Matrix.py

#Clean up
#rm $inputsPath"/Format_Matrix.py"

#Print status message
echo "Analysis complete!"

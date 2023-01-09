#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_SNP_jobOutput
#$ -pe smp 8
#Script to run the SSR pipeline
#Usage: qsub ssr_pipeline_SNP.sh inputsFile
#Usage Ex: qsub ssr_pipeline_SNP.sh inputPaths_romero_run1.txt
#Usage Ex: qsub ssr_pipeline_SNP.sh inputPaths_romero_run2.txt
#Usage Ex: qsub ssr_pipeline_SNP.sh inputPaths_romero_run3.txt
#Usage Ex: qsub ssr_pipeline_SNP.sh inputPaths_romero_run4.txt
#Usage Ex: qsub ssr_pipeline_SNP.sh inputPaths_romero_run5.txt

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

# setup the inputs path
inputsPath=$outputsPath"/"$projectDir"_SSR_prep"

#Make a new directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_SNP"
mkdir $outputsPath
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile=$outputsPath"/pipeline_SNP_summary.txt"
#Add pipeline info to outputs
echo -e "SSR pipline SNP inputs for $projectDir \n" > $inputOutFile


#SSR Analysis Stage - SNP Calling Workflow

#Set input paths
inputsPath=$inputsPath"/aligned"

#Copy pipeline scripts to inputs directory
cp SamIAm.py $inputsPath
cp Format_VCF-Matrix.py $inputsPath

#Move to the inputs directory
cd $inputsPath

#Loop through all filtered sam files
for f1 in "$inputsPath"/*.sam; do
	#Print status message
	echo "Processing $f1"
	#Run SSR pipeline
	python2 SamIAm.py -sam $f1 -p "yes"
	#Write inputs out to summary file
	echo "python2 SamIAm.py -sam $f1 -p yes" >> $inputOutFile
	#Status message
	echo "Processed!"
done

#Retrieve and format sample tag list
sampleTags=$(for i in "$inputsPath"/*.sam; do basename $i | sed "s/^/\"/g" | sed "s/\.sam/\",/g" | tr '\n' ' '; done)
sampleTags=$(echo $sampleTags | sed 's/.$//')

#Find and replace the sample list
sed -i "s/\"FIND_ME_REPLACE_ME\"/$sampleTags/g" Format_VCF-Matrix.py

#Format matrix
python2 Format_VCF-Matrix.py
echo "python2 Format_VCF-Matrix.py" >> $inputOutFile

#Clean up
rm $inputsPath"/SamIAm.py"
rm $inputsPath"/Format_VCF-Matrix.py"

#Re-name and move output matrix
mv VCF_Matrix.txt $outputsPath"/"$projectDir"_VCF_Matrix.txt"

#Print status message
echo "Analysis complete!"

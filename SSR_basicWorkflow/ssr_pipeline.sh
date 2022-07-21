#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_pipeline_jobOutput
#$ -pe smp 8
#Script to run the SSR pipeline
#Usage: qsub ssr_pipeline.sh inputsFile
#Usage Ex: qsub ssr_pipeline.sh inputPaths_romero_run1.txt
#Usage Ex: qsub ssr_pipeline.sh inputPaths_romero_run2.txt
#Usage Ex: qsub ssr_pipeline.sh inputPaths_romero_run3.txt
#Usage Ex: qsub ssr_pipeline.sh inputPaths_romero_run4.txt
#Usage Ex: qsub ssr_pipeline.sh inputPaths_romero_run5.txt

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
readPath=$(grep "pairedReads:" ../"InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")
#Retrieve adapter absolute path for alignment
infoPath=$(grep "info:" ../"InputData/"$inputsFile | tr -d " " | sed "s/info://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

#Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir"_SSR_basicWorkflow"
mkdir $outputsPath
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile=$outputsPath"/pipeline_summary.txt"
versionFile=$outputsPath"/software_summary.txt"
#Add pipeline info to outputs
echo -e "SSR pipline inputs for $projectDir \n" > $inputOutFile
echo -e "SSR pipeline software versions for $projectDir \n" > $versionFile


#Read Prep Stage

#Move to directory with analysis prep scripts
cd ../Prep
#Quality control with fastqc
bash fastqc_ssr_projects.sh $1
#Trimming with trimmomatic
bash trimmomatic_ssr_projects.sh $1
#Mapping with bwa
bash bwa_ssr_projects.sh $1


#SSR Analysis Stage

#Set input paths
inputsPath=$outputsPath"/aligned"

#Move to directory with SSR analysis scripts
cd ../SSR_basicWorkflow
#Copy pipeline scripts to inputs directory
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
	#Write inputs out to summary file
	echo "python2 GapGenes.v3.py -sam $f1 -C $infoPath -P unpaired" >> $inputOutFile
	echo "python2 SnipMatrix.py "$f1".Matrix.txt" >> $inputOutFile
	#Status message
	echo "Processed!"
done

#Retrieve and format sample tag list
sampleTags=$(for i in "$inputsPath"/*.sam; do basename $i | sed "s/^/\"/g" | sed "s/\.sam/\",/g" | tr '\n' ' ' | sed 's/..$//'; done)

#Find and replace the sample list
sed -i "s/\"FIND_ME_REPLACE_ME\"/$sampleTags/g" "$inputsPath"/Format_Matrix.py

#Format matrix
python2 Format_Matrix.py
echo "python2 Format_Matrix.py" >> $inputOutFile

#Clean up
rm $inputsPath"/GapGenes.v3.py"
rm $inputsPath"/SnipMatrix.py"
rm $inputsPath"/Format_Matrix.py"

#Print status message
echo "Analysis complete!"

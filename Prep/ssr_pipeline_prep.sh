#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_prep_jobOutput
#$ -pe smp 8
#Script to run the SSR pipeline
#Usage: qsub ssr_pipeline_prep.sh inputsFile
#Usage Ex: qsub ssr_pipeline_prep.sh inputPaths_romero_run1.txt
#Usage Ex: qsub ssr_pipeline_prep.sh inputPaths_romero_run2.txt
#Usage Ex: qsub ssr_pipeline_prep.sh inputPaths_romero_run3.txt
#Usage Ex: qsub ssr_pipeline_prep.sh inputPaths_romero_run4.txt
#Usage Ex: qsub ssr_pipeline_prep.sh inputPaths_romero_run5.txt

#Required modules for ND CRC servers
module load bio

#Retrieve input argument of a inputs file
inputsFile=$1

#Retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
#Retrieve genome reference absolute path for alignment
ref=$(grep "genomeReference:" ../"InputData/"$inputsFile | tr -d " " | sed "s/genomeReference://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

#Make a new directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_prep"
mkdir $outputsPath
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
versionFile=$outputsPath"/software_prep_summary.txt"
#Add pipeline info to outputs
echo -e "SSR pipeline prep software versions for $projectDir \n" > $versionFile


#Analysis Prep Stage

#Make sure the reference genome has been indexed
bwa index $ref

#Quality control with fastqc
bash fastqc_ssr_projects.sh $1
#Trimming with trimmomatic
bash trimmomatic_ssr_projects.sh $1
#Mapping with bwa
bash bwa_ssr_projects.sh $1

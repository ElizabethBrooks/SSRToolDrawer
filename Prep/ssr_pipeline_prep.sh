#!/bin/bash

#Script to run the SSR pipeline
#Usage: bash ssr_pipeline_prep.sh inputsFile outputsPath

#Required modules for ND CRC servers
module load bio

#Retrieve input argument of a inputs file
inputsFile=$1

# retrieve input outputs path
outputsPath=$2

#Retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
#Retrieve genome reference absolute path for alignment
ref=$(grep "genomeReference:" ../"InputData/"$inputsFile | tr -d " " | sed "s/genomeReference://g")

#Name output file of inputs
versionFile=$outputsPath"/software_prep_summary.txt"
#Add pipeline info to outputs
echo -e "SSR pipeline prep software versions for $projectDir \n" > $versionFile


#Analysis Prep Stage

#Make sure the reference genome has been indexed
bwa index $ref

#Quality control with fastqc
bash fastqc_ssr_projects.sh $inputsFile
#Trimming with trimmomatic
bash trimmomatic_ssr_projects.sh $inputsFile
#Mapping with bwa
bash bwa_ssr_projects.sh $inputsFile

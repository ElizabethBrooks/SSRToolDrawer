#!/bin/bash

# script to prep input data for the SSR pipelines
# usage: bash ssr_pipeline_prep.sh inputsFile outputsPath

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve analysis outputs absolute path
outputsPath=$2

# retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
# retrieve reference absolute path for alignment
ref=$(grep "reference:" ../"InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")

# name output file of inputs
versionFile=$outputsPath"/software_prep_summary.txt"
# add pipeline info to outputs
echo -e "SSR pipeline prep software versions for $projectDir \n" > $versionFile


# Analysis Prep Stage

# status message
echo "Prep started..."

# make sure the ssr info file has been converted from csv to txt
#bash ssr_info_prep.sh

# make sure the reference has been indexed
#bwa index $ref

# move to scripts directory
cd Scripts

# quality control with fastqc
bash qc_fastqc.sh $inputsFile $outputsPath
# trimming with trimmomatic
bash trimming_trimmomatic.sh $inputsFile $outputsPath
# mapping with bwa
bash alignment_bwa.sh $inputsFile $outputsPath

# status message
echo "Prep complete!"

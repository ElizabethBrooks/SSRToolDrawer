#!/bin/bash

# script to run the SSR pipeline
# usage: bash ssr_pipeline_prep.sh runInputs outputsPath

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve input outputs path
outputsPath=$2

# retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
# retrieve genome reference absolute path for alignment
ref=$(grep "genomeReference:" ../"InputData/inputPaths_ssr_pipeline.txt" | tr -d " " | sed "s/genomeReference://g")

# name output file of inputs
versionFile=$outputsPath"/software_prep_summary.txt"
# add pipeline info to outputs
echo -e "SSR pipeline prep software versions for $projectDir \n" > $versionFile


# Analysis Prep Stage

# status message
echo "Prep started..."

# make sure the ssr info file has been converted from csv to txt
#bash ssr_info_prep.sh

# make sure the reference genome has been indexed
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

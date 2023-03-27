#!/bin/bash

# script to prep input data for the SSR pipelines
# usage: bash ssr_pipeline_prep.sh inputsPath inputsFile baseDir

# retrieve input argument of a inputs file
inputsPath=$1

# retrieve analysis outputs absolute path
inputsFile=$2

# retrieve base of working directory
baseDir=$3

# retrieve the project ID 
projectDir=$(grep "ID:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/ID://g")

# name output file of inputs
versionFile=$inputsPath"/software_prep_summary.txt"
# add pipeline info to outputs
echo -e "SSR pipeline prep software versions for $projectDir \n" > $versionFile


# Analysis Prep Stage

# status message
echo "Prep started..."

# move to scripts directory
cd Scripts

# quality control with fastqc
#bash qc_fastqc.sh $inputsFile $inputsPath $baseDir
# trimming with trimmomatic
bash trimming_trimmomatic.sh $inputsFile $inputsPath $baseDir
# mapping with bwa
bash alignment_bwa.sh $inputsFile $inputsPath $baseDir

# status message
echo "Prep complete!"

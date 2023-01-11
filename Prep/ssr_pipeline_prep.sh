#!/bin/bash

# script to run the SSR pipeline
# usage: bash ssr_pipeline_prep.sh inputsFile inputsPath

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve input outputs path
inputsPath=$2

# retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
# retrieve genome reference absolute path for alignment
ref=$(grep "genomeReference:" ../"InputData/"$inputsFile | tr -d " " | sed "s/genomeReference://g")

# name output file of inputs
versionFile=$inputsPath"/software_prep_summary.txt"
# add pipeline info to outputs
echo -e "SSR pipeline prep software versions for $projectDir \n" > $versionFile


# Analysis Prep Stage

# status message
echo "Prep started..."

# make sure the reference genome has been indexed
#bwa index $ref

# quality control with fastqc
bash qc_fastqc.sh $inputsFile $inputsPath
# trimming with trimmomatic
bash trimming_trimmomatic.sh $inputsFile $inputsPath
# mapping with bwa
bash alignment_bwa.sh $inputsFile $inputsPath

# status message
echo "Prep complete!"

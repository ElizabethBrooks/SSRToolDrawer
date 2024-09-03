#!/bin/bash

# script to prep input data for the SSR pipelines
# usage: bash ssr_pipeline_prep.sh inputsPath inputsFile baseDir

# retrieve input argument of a inputs file
inputsPath=$1

# retrieve analysis outputs absolute path
inputsFile=$2

# retrieve base of working directory
baseDir=$3

# retrieve analysis outputs path
outputsPath=$(grep "outputs:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/outputs://g")

# retrieve analysis outputs path
outputsPath=$(dirname $inputsPath)

# name output file of inputs
versionFile=$outputsPath"/info/software_summary_prep.txt"

# add pipeline info to outputs
echo -e "SSR pipeline prep software versions: \n" > $versionFile

# setup reports outputs path
reportsOut=$outputsPath"/reports"

# make reports directory
mkdir $reportsOut


# Analysis Prep Stage

# status message
echo "Prep started..."

# move to scripts directory
cd Scripts

# quality control with fastqc
bash qc_fastqc.sh $inputsPath $inputsFile $baseDir

# create qc report with multiqc
bash report_multiqc.sh $inputsPath $inputsFile $baseDir "raw"

# trimming with trimmomatic
bash trimming_trimmomatic.sh $inputsPath $inputsFile $baseDir

# create qc report with multiqc
bash report_multiqc.sh $inputsPath $inputsFile $baseDir "trimmed"

# mapping with bwa
bash alignment_bwa.sh $inputsPath $inputsFile $baseDir

# status message
echo "Prep complete!"

#!/bin/bash

# script to run the SSR pipeline
# usage: bash ssr_pipeline_prep.sh runInputs outputsPath
# usage Ex: qsub ssr_pipeline_basic.sh inputPaths_run1.txt /afs/crc.nd.edu/group/genomics/Mando/GBCF_bioinformatics_romero_SSR
# usage Ex: qsub ssr_pipeline_basic.sh inputPaths_run2.txt /afs/crc.nd.edu/group/genomics/Mando/GBCF_bioinformatics_romero_SSR
# usage Ex: qsub ssr_pipeline_basic.sh inputPaths_run3.txt /afs/crc.nd.edu/group/genomics/Mando/GBCF_bioinformatics_romero_SSR
# usage Ex: qsub ssr_pipeline_basic.sh inputPaths_run4.txt /afs/crc.nd.edu/group/genomics/Mando/GBCF_bioinformatics_romero_SSR
# usage Ex: qsub ssr_pipeline_basic.sh inputPaths_run5.txt /afs/crc.nd.edu/group/genomics/Mando/GBCF_bioinformatics_romero_SSR
# usage Ex: qsub ssr_pipeline_basic.sh inputPaths_run6.txt /afs/crc.nd.edu/group/genomics/Mando/GBCF_bioinformatics_romero_SSR

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve input outputs path
outputsPath=$2

# retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
# retrieve reference absolute path for alignment
ref=$(grep "reference:" ../"InputData/inputPaths_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")

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

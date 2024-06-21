#!/bin/bash

# script to perform trimming of paired end reads
# usage: bash report_multiqc.sh inputsPath inputType

# required modules for ND CRC servers
#module load bio

# retrieve input argument of a inputs file
inputsPath=$1

# retrieve inputs type
inputType=$2

# retrieve analysis outputs path
outputsPath=$(dirname $inputsPath)

# name output file of inputs
versionFile=$outputsPath"/info/software_summary_prep.txt"

# output software version
echo "QC report: " >> $versionFile
multiqc --version >> $versionFile
echo -e "\n" >> $versionFile

# setup multiqc outputs path
multiqcOut=$outputsPath"/reports/qc_"$inputType

# make multiqc outputs directory
mkdir $multiqcOut

# name of directory with qc results
qcOut=$inputsPath"/qc"

# move to the new directory
cd $qcOut

# run multiqc
multiqc $qcOut

# move multiqc outputs
mv multiqc* $multiqcOut

# print status message
echo "Analysis complete!"

#!/bin/bash

# script to perform trimming of paired end reads
# usage: bash report_multiqc.sh inputsPath inputsFile baseDir inputType

# required modules for ND CRC servers
#module load bio

# retrieve input argument of a inputs file
inputsPath=$1

# retrieve analysis outputs absolute path
inputsFile=$2

# retrieve base of working directory
baseDir=$3

# retrieve inputs type
inputType=$4

# retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/pairedReads://g")

# retrieve analysis outputs path
outputsPath=$(dirname $inputsPath)

# name output file of inputs
versionFile=$outputsPath"/info/software_summary_prep.txt"

# output software version
echo "QC report: " >> $versionFile
multiqc --version >> $versionFile
echo -e "\n" >> $versionFile

# retrieve the run number 
runNum=$(grep "run:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/run://g")

# setup multiqc outputs path
multiqcOut=$outputsPath"/reports/qc_"$inputType"_"$runNum

# this is used to re-start analyses
# check if the outputs directory already exsists
if [ -d $multiqcOut ]; then
	echo "$multiqcOut directory exists... exiting"
	exit
fi

# make the new directory
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

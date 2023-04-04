#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N test_ssr_SNP_jobOutput
#$ -pe smp 4

# script to create a formatted variant matrix
# usage: bash variant_matrix_formatting.sh inputsPath inputsFile baseDir


# retrieve input outputs path
inputsPath=$1

# retrieve inputs file
inputsFile=$2

# retrieve base of working directory
baseDir=$3

# retrieve the run number 
runNum=$(grep "run:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/run://g")
# retrieve the project ID 
projectDir=$(grep "ID:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
# retrieve analysis outputs path
outputsPath=$(grep "outputs:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/outputs://g")
# retrieve ssr regions path
regionsPath=$(grep "ssrRegions:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrRegions://g")

# set outputs path
outputsPath=$outputsPath"/"$projectDir"_SSR_SNP"

# set the inputs path
inputsPath=$inputsPath"/variants"

# set output matrix file
resultsFile=$outputsPath"/"$runNum".txt"

# contig (marker) ID list
contingList=$(cat $regionsPath | cut -f 1)

# loop over each sample
for f1 in $inputsPath"/variants/"*".noHeader.vcf"; do
	# retrieve sample tag
	sampleTag=$(basename $f1 | sed "s/\.noHeader\.vcf$//g")
	# retrieve 
done



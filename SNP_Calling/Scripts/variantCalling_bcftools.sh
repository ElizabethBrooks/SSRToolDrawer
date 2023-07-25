#!/bin/bash

# script to run the SSR pipeline
# usage: bash variantCalling_bcftools.sh inputsFile baseDir inputsPathList

# retrieve inputs file
inputsFile=$2

# retrieve base of working directory
baseDir=$3

# set run num tag
runNum="combined"
# retrieve reference path
ref=$(grep "reference:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")
# retrieve ssr regions path
regionsPath=$(grep "ssrRegions:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrRegions://g")
# retrieve primers path
primerPath=$(grep "primers:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/primers://g")

# setup the variant calling directory
outputsPath=$inputsPath"/variantsCalled"
# create the directory
mkdir $outputsPath

# name of output file of inputs
versionFile=$inputsPath"/software_VC_summary.txt"

# output software version
echo "Variant calling: " >> $versionFile
bcftools --version >> $versionFile

# set inputs inputsPath
inputsPath=$inputsPath"/clipped"


# Variant Calling Stage - SNP Calling Workflow

# pre-clean up
rm $outputsPath"/inputBAMList.txt"

# loop over each input file
for arg in "${@:3}"; do
	# retrieve input outputs path
	inputsPath=$1
	# loop through all filtered sam files
	ls $inputsPath"/"*".readGroups.bam" >> $outputsPath"/inputBAMList.txt"
done

# status message
echo "Performing variant calling..."

# TO-DO
# https://github.com/samtools/bcftools/issues/811
# bcftools mpileup -a AD
# bcftools call -G

# calculate the read coverage of positions in the genome
bcftools mpileup --threads 4 -d 8000 -f $ref -Ob -o $outputsPath"/"$runNum"_raw.bcf" -b $outputsPath"/inputBAMList.txt"
rm $outputsPath"/inputBAMList.txt"
# detect the single nucleotide polymorphisms 
bcftools call --threads 4 -mv -Oz -o $outputsPath"/"$runNum"_calls.vcf.gz" $outputsPath"/"$runNum"_raw.bcf" 
rm $outputsPath"/"$runNum"_raw.bcf"
# index vcf file
bcftools index --threads 4 $outputsPath"/"$runNum"_calls.vcf.gz"
# normalize indels
bcftools norm --threads 4 -f $ref -Ob -o $outputsPath"/"$runNum"_calls.norm.bcf" $outputsPath"/"$runNum"_calls.vcf.gz"
rm $outputsPath"/"$runNum"_calls.vcf.gz"*

# status message
echo "Analysis conplete!"

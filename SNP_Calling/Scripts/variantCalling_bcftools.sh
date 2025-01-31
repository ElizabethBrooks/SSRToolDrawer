#!/bin/bash

# script to run the SSR pipeline
# usage: bash variantCalling_bcftools.sh baseDir inputsPath

# retrieve base of working directory
baseDir=$1

# set run num tag
runNum="combined"

# retrieve analysis inputs path
inputsPath=$2

# retrieve reference path
ref=$(grep "reference:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")

# retrieve ssr regions path
regionsPath=$(grep "ssrRegions:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrRegions://g")

# retrieve primers path
primerPath=$(grep "primers:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/primers://g")

# retrieve analysis outputs path
outputsPath=$(dirname $inputsPath)

# name output file of inputs
versionFile=$outputsPath"/info/software_summary_SNP.txt"

# output software version
echo "Variant calling: " >> $versionFile
bcftools --version >> $versionFile
echo -e "\n" >> $versionFile

# setup the variant calling directory
outputsPath=$inputsPath"/variantsCalled"
# create the directory
mkdir $outputsPath


# Variant Calling Stage - SNP Calling Workflow

# add the inputs file paths to a txt file
ls -d $inputsPath"/SNP_Calling_prep_run"*"/clipped/"*".readGroups.bam" | grep -v "Undetermined" > $outputsPath"/inputBAMList.txt"

# status message
echo "Performing variant calling..."

## TO-DO
## https://github.com/samtools/bcftools/issues/811
## bcftools mpileup -a AD
## bcftools call -G

# calculate the read coverage of positions in the genome
bcftools mpileup --threads 8 -d 8000 -f $ref -Ob -o $outputsPath"/"$runNum"_raw.bcf" -b $outputsPath"/inputBAMList.txt"
#rm $outputsPath"/inputBAMList.txt"
# detect the single nucleotide polymorphisms 
bcftools call --threads 8 -mv -Oz -o $outputsPath"/"$runNum"_calls.vcf.gz" $outputsPath"/"$runNum"_raw.bcf" 
rm $outputsPath"/"$runNum"_raw.bcf"
# index vcf file
bcftools index --threads 8 $outputsPath"/"$runNum"_calls.vcf.gz"
# normalize indels
bcftools norm --threads 8 -f $ref -Ob -o $outputsPath"/"$runNum"_calls.norm.bcf" $outputsPath"/"$runNum"_calls.vcf.gz"
rm $outputsPath"/"$runNum"_calls.vcf.gz"*

# status message
echo "Analysis conplete!"

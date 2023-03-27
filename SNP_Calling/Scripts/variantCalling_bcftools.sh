#!/bin/bash

# script to run the SSR pipeline
# usage: bash variantCalling_bcftools.sh inputsPath baseDir inputsFile

# retrieve input outputs path
inputsPath=$1

# retrieve base of working directory
baseDir=$2

# retrieve inputs file
inputsFile=$3

# retrieve the run number 
runNum=$(grep "run:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/run://g")
# retrieve the project ID 
projectDir=$(grep "ID:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
# retrieve reference path
ref=$(grep "reference:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")
# retrieve ssr regions path
regionsPath=$(grep "ssrRegions:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrRegions://g")

# setup the variant calling directory
outputsPath=$inputsPath"/variants"
# create the directory
mkdir $outputsPath

# name of output file of inputs
versionFile=$inputsPath"/software_VC_summary.txt"

# retrieve software version
bcftools --version >> $versionFile


# Variant Calling Stage - SNP Calling Workflow

# loop through all filtered sam files
ls $inputsPath"/clipped/"*".readGroups.bam" > $outputsPath"/inputBAMList.txt"

# status message
echo "Performing variant calling for $runNum"

# calculate the read coverage of positions in the genome
bcftools mpileup --threads 4 -d 8000 -f $ref -Ob -o $outputsPath"/"$runNum"_raw.bcf" -b $outputsPath"/inputBAMList.txt"
rm $outputsPath"/inputBAMList.txt"
# detect the single nucleotide polymorphisms 
bcftools call --threads 4 -mv -Oz -o $outputsPath"/"$runNum"_calls.vcf.gz" $outputsPath"/"$runNum"_raw.bcf" 
rm $outputsPath"/"$runNum"_raw.bcf"
# index vcf file
bcftools index --threads 4 $outputsPath"/"$runNum"_calls.vcf.gz"
# normalize indels
bcftools norm --threads 4 -f $ref -o $outputsPath"/"$runNum"_calls.norm.bcf" $outputsPath"/"$runNum"_calls.vcf.gz"
rm $outputsPath"/"$runNum"_calls.vcf.gz"*
# convert from BCF to VCF
bcftools view --threads 4 -Ov -o $outputsPath"/"$runNum"_calls.norm.vcf" $outputsPath"/"$runNum"_calls.norm.bcf"
rm $outputsPath"/"$runNum"_calls.norm.vcf"
# remove ssr regions
bedtools intersect -v -a $outputsPath"/"$runNum"_calls.norm.vcf" -b $regionsPath -header > $outputsPath"/"$runNum"_noSSR.vcf"
# TO-DO
# remove primer regions

# status message
echo "Analysis conplete!"

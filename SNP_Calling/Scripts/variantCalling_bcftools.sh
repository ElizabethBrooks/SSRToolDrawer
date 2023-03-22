#!/bin/bash

# script to run the SSR pipeline
# usage: bash variantCalling_bcftools.sh inputsPath projectDir reference runNum

# retrieve input outputs path
inputsPath=$1

# retrieve the project ID 
projectDir=$2

# retrieve genome reference absolute path for alignment
ref=$3

# retrieve the run number
runNum=$4

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
#rm $outputsPath"/inputBAMList.txt"
# detect the single nucleotide polymorphisms 
bcftools call --threads 4 -mv -Oz -o $outputsPath"/"$runNum"_calls.vcf.gz" $outputsPath"/"$runNum"_raw.bcf" 
#rm $outputsPath"/"$runNum"_raw.bcf"
# index vcf file
bcftools index --threads 4 $outputsPath"/"$runNum"_calls.vcf.gz"
# normalize indels
bcftools norm --threads 4 -f $ref -o $outputsPath"/"$runNum"_calls.norm.bcf" $outputsPath"/"$runNum"_calls.vcf.gz"
#rm $outputsPath"/"$runNum"_calls.vcf.gz"*
# convert from BCF to VCF
bcftools view --threads 4 -Ov -o $outputsPath"/"$runNum"_calls.norm.vcf" $outputsPath"/"$runNum"_calls.norm.bcf"
#rm $outputsPath"/"$runNum"_calls.norm.vcf"

# status message
echo "Analysis conplete!"

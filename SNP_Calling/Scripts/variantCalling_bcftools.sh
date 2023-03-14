#!/bin/bash

# script to run the SSR pipeline
# usage: bash variantCalling_bcftools.sh inputsFile inputsPath

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve input outputs path
inputsPath=$2

# retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
# retrieve genome reference absolute path for alignment
ref=$(grep "genomeReference:" ../"InputData/"$inputsFile | tr -d " " | sed "s/genomeReference://g")

# setup the variant calling directory
dataPath=$inputsPath"/variants"
# create the directory
mkdir $dataPath

# name of output file of inputs
versionFile=$inputsPath"/software_VC_summary.txt"

# retrieve software version
bcftools --version >> $versionFile


# Variant Calling Stage - SNP Calling Workflow

# set input paths
inputsPath=$inputsPath"/sorted"

# loop through all filtered sam files
ls $inputsPath"/"*"sortedCoordinate.bam" > "inputBAMList.txt"

# remove the file path
noPath=$(echo $inputsFile | sed 's/\.txt$//g' | sed 's/inputPaths_//g')

# status message
echo "Performing variant calling for "$noPath"..."

# calculate the read coverage of positions in the genome
bcftools mpileup --threads 4 -d 8000 -f $ref -Ob -o $dataPath"/"$noPath"_raw.bcf" -b "inputBAMList.txt"
# detect the single nucleotide polymorphisms 
bcftools call --threads 4 -mv -Oz -o $dataPath"/"$noPath"_calls.vcf.gz" $dataPath"/"$noPath"_raw.bcf" 
# index vcf file
bcftools index --threads 4 $dataPath"/"$noPath"_calls.vcf.gz"
# normalize indels
bcftools norm --threads 4 -f $ref -o $dataPath"/"$noPath"_calls.norm.bcf" $dataPath"/"$noPath"_calls.vcf.gz"
# filter adjacent indels within 5bp
bcftools filter --threads 4 --IndelGap 5 -Ob -o $dataPath"/"$noPath"_calls.norm.flt-indels.bcf" $dataPath"/"$noPath"_calls.norm.bcf"
# convert from BCF to VCF
bcftools view --threads 4 -Ov -o $dataPath"/"$noPath"_calls.norm.flt-indels.vcf" $dataPath"/"$noPath"_calls.norm.flt-indels.bcf"

# status message
echo "Analysis conplete!"

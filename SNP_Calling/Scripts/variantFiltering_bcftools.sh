#!/bin/bash

# script to run the SSR pipeline
# usage: bash variantFiltering_bcftools.sh baseDir

# retrieve base of working directory
baseDir=$1

# set the run number 
runNum="combined"

# retrieve analysis inputs path
inputsPath=$(grep "outputs:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/outputs://g")
inputsPath=$inputsPath"/SSR_SNP"
# retrieve reference path
ref=$(grep "reference:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")

# setup the variant calling directory
outputsPath=$inputsPath"/variantsFiltered"
# create the directory
mkdir $outputsPath

# name of output file of inputs
versionFile=$outputsPath"/software_VF_summary.txt"

# name of output file of filtering summary
outputsFile=$outputsPath"/variantFiltering_stats.txt"

# set inputs path
inputsPath=$inputsPath"/variantsCalled"

# output software version
echo "Variant filtering: " >> $versionFile
bcftools --version >> $versionFile


# Variant Filtering Stage - SNP Calling Workflow

# status message
echo "Performing variant filtering for $runNum"

#Check total variants
echo "Total variants from filtered reads with MQ > 60: " > $outputsFile
bcftools view $inputsPath"/"$runNum"_calls.norm.bcf" | grep -v "#" | wc -l >> $outputsFile

#Include sites with quality > 20 
bcftools filter --threads 4 -i '%QUAL>20' $inputsPath"/"$runNum"_calls.norm.bcf" -Ob -o $outputsPath"/"$runNum"_calls.flt-qual.bcf"
echo "& including sites with quality > 20: " >> $outputsFile
bcftools view --threads 4 $outputsPath"/"$runNum"_calls.flt-qual.bcf" | grep -v "#" | wc -l >> $outputsFile

#Include sites with average read depth > 10
bcftools filter --threads 4 -i 'INFO/DP>10' $outputsPath"/"$runNum"_calls.flt-qual.bcf" -Ob -o $outputsPath"/"$runNum"_calls.flt-qualDP.bcf"
echo "& including sites with average read depth > 10: " >> $outputsFile
bcftools view --threads 4 $outputsPath"/"$runNum"_calls.flt-qualDP.bcf" | grep -v "#" | wc -l >> $outputsFile

#Turn on left alignment, normalize indels
## and split multiallelic sites
##bcftools norm --threads 4 -m -any -f $ref $outputsPath"/"$runNum"_calls.flt-qualDP.bcf" -Ov -o $outputsPath"/"$runNum"_calls.flt-norm.vcf"
bcftools norm --threads 4 -f $ref $outputsPath"/"$runNum"_calls.flt-qualDP.bcf" -Ov -o $outputsPath"/"$runNum"_calls.flt-norm.vcf"
echo "& with left alignment and normalized indels: " >> $outputsFile
bcftools view --threads 4 $outputsPath"/"$runNum"_calls.flt-norm.vcf" | grep -v "#" | wc -l >> $outputsFile

# index the bcf
#bcftools index $outputsPath"/"$runNum"_calls.flt-norm.vcf"

# clean up
#rm $outputsPath"/"$runNum"_calls.flt-qual.bcf"
#rm $outputsPath"/"$runNum"_calls.flt-qualDP.bcf"

# status message
echo "Analysis conplete!"

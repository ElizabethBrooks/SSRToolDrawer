#!/bin/bash

# script to run the SSR pipeline
# usage: bash variantTrimming_bedtools.sh baseDir inputsPath

# retrieve base of working directory
baseDir=$1

# set run num tag
runNum="combined"

# retrieve analysis inputs path
inputsPath=$2

# retrieve ssr regions path
regionsPath=$(grep "ssrRegions:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrRegions://g")

# retrieve primers path
primerPath=$(grep "primers:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/primers://g")

# retrieve analysis outputs path
outputsPath=$(dirname $inputsPath)

# name output file of inputs
versionFile=$outputsPath"/info/software_summary_SNP.txt"

# setup the variant calling directory
outputsPath=$inputsPath"/variantsTrimmed"
# create the directory
mkdir $outputsPath

# set inputs path
inputsPath=$inputsPath"/variantsFiltered"

# output software version
echo "Variant trimming: " >> $versionFile
bcftools --version >> $versionFile
echo -e "\n" >> $versionFile


# Variant Trimming Stage - SNP Trimming Workflow

# status message
echo "Performing variant trimming for $runNum"

# remove ssr regions
bedtools intersect -v -header -a $inputsPath"/"$runNum"_calls.flt-norm.vcf" -b $regionsPath > $outputsPath"/"$runNum"_noSSR.vcf"
#rm $outputsPath"/"$runNum"_calls.norm.vcf"
# remove sense primer regions
cat $primerPath | cut -f 1-3 > $outputsPath"/tmp_sense_primerRegions.bed"
bedtools intersect -v -header -a $outputsPath"/"$runNum"_noSSR.vcf" -b $outputsPath"/tmp_sense_primerRegions.bed" > $outputsPath"/"$runNum"_noSense.vcf"
rm $outputsPath"/tmp_sense_primerRegions.bed"
rm $outputsPath"/"$runNum"_noSSR.vcf"
# remove antisense primer regions
cat $primerPath | cut -f 4-6 > $outputsPath"/tmp_antisense_primerRegions.bed"
bedtools intersect -v -header -a $outputsPath"/"$runNum"_noSense.vcf" -b $outputsPath"/tmp_antisense_primerRegions.bed" > $outputsPath"/"$runNum"_trimmed.vcf"
rm $outputsPath"/tmp_antisense_primerRegions.bed"
rm $outputsPath"/"$runNum"_noSense.vcf"

# status message
echo "Analysis conplete!"

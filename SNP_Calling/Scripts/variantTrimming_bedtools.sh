#!/bin/bash

# script to run the SSR pipeline
# usage: bash variantTrimming_bedtools.sh inputsPath inputsFile baseDir

# retrieve input outputs path
inputsPath=$1

# retrieve inputs file
inputsFile=$2

# retrieve base of working directory
baseDir=$3

# retrieve the run number 
runNum=$(grep "run:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/run://g")
# retrieve ssr regions path
regionsPath=$(grep "ssrRegions:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrRegions://g")
# retrieve primers path
primerPath=$(grep "primers:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/primers://g")

# setup the variant calling directory
outputsPath=$inputsPath"/variants"

# name of output file of inputs
versionFile=$inputsPath"/software_VC_summary.txt"

# retrieve software version
bcftools --version >> $versionFile


# Variant Trimming Stage - SNP Trimming Workflow

# status message
echo "Performing variant trimming for $runNum"

# remove ssr regions
bedtools intersect -v -header -a $outputsPath"/"$runNum"_calls.norm.vcf" -b $regionsPath > $outputsPath"/"$runNum"_noSSR.vcf"
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

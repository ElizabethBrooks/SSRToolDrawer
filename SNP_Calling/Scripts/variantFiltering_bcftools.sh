#!/bin/bash

# script to run the SSR pipeline
# usage: bash variantFiltering_bcftools.sh inputsPath inputsFile baseDir

# retrieve input outputs path
inputsPath=$1

# retrieve inputs file
inputsFile=$2

# retrieve base of working directory
baseDir=$3

# retrieve the run number 
runNum=$(grep "run:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/run://g")
# retrieve reference path
ref=$(grep "reference:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")

# setup the variant calling directory
outputsPath=$inputsPath"/variantsFiltered"
# create the directory
mkdir $outputsPath

# name of output file of inputs
versionFile=$inputsPath"/software_VC_summary.txt"

# name of output file of filtering summary
outputsFile=$inputsPath"/filtering_summary.txt"

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
bcftools filter --threads 4 -i '%QUAL<1001' $inputsPath"/"$runNum"_calls.norm.bcf" | grep -v "#" | wc -l >> $outputsFile

#Include sites with quality > 20 
bcftools filter --threads 4 -i '%QUAL>20' $inputsPath"/"$runNum"_calls.norm.bcf" -Ob -o $outputsPath"/"$runNum"_calls.flt-qual.bcf"
echo "bcftools filter --threads 4 -i '%QUAL>20' "$inputsPath"/"$runNum"_calls.norm.bcf -Ob -o "$outputsPath"/"$runNum"_calls.flt-qual.bcf" >> $inputOutFile
echo "& including sites with quality > 20: " >> $outputsFile
bcftools filter --threads 4 -i '%QUAL>20' $inputsPath"/"$runNum"_calls.norm.bcf" | grep -v "#" | wc -l >> $outputsFile

#Include sites with average read depth > 10
bcftools filter --threads 4 -i 'INFO/DP>10' $outputsPath"/"$runNum"_calls.flt-qual.bcf" -Ob -o $outputsPath"/"$runNum"_calls.flt-qualDP.bcf"
echo "bcftools filter --threads 4 -i 'INFO/DP>10' "$outputsPath"/"$runNum"_calls.flt-qual.bcf -Ob -o "$outputsPath"/"$runNum"_calls.flt-qualDP.bcf" >> $inputOutFile
echo "& including sites with average read depth > 10: " >> $outputsFile
bcftools filter --threads 4 -i 'INFO/DP>10' $outputsPath"/"$runNum"_calls.flt-qual.bcf" | grep -v "#" | wc -l >> $outputsFile

#Exclude heterozygous sites
bcftools filter --threads 4 -e 'GT="het"' $outputsPath"/"$runNum"_calls.flt-qualDP.bcf" -Ob -o $outputsPath"/"$runNum"_calls.flt-qualDP-homo.bcf"
echo "bcftools filter --threads 4 -e 'GT=\"het\"' "$outputsPath"/"$runNum"_calls.flt-qualDP.bcf -Ob -o "$outputsPath"/"$runNum"_calls.flt-qualDP-homo.bcf" >> $inputOutFile
echo "& excluding heterozygous sites: " >> $outputsFile
bcftools filter --threads 4 -e 'GT="het"' $outputsPath"/"$runNum"_calls.flt-qualDP.bcf" | grep -v "#" | wc -l >> $outputsFile

#Turn on left alignment, normalize indels, and collapse multi allelic sites
bcftools norm --threads 4 -m +any -f $ref $outputsPath"/"$runNum"_calls.flt-qualDP-homo.bcf" -Ob -o $outputsPath"/"$runNum"_calls.flt-norm.bcf"
echo "bcftools norm --threads 4 -m +any -f "$ref" "$outputsPath"/"$runNum"_calls.flt-qualDP-homo.bcf -Ob -o "$outputsPath"/"$runNum"_calls.flt-norm.bcf" >> $inputOutFile
echo "& with left alignment, normalized indels, and collapsed multi allelic sites: " >> $outputsFile
bcftools norm --threads 4 -m +any -f $ref $outputsPath"/"$runNum"_calls.flt-qualDP-homo.bcf" | grep -v "#" | wc -l >> $outputsFile

#Index bcf file
bcftools index --threads 4 $outputsPath"/"$runNum"_calls.flt-norm.bcf"
echo "bcftools index --threads 4 "$outputsPath"/"$runNum"_calls.flt-norm.bcf" >> $inputOutFile

# clean up
#rm $outputsPath"/"$runNum"_calls.flt-qual.bcf"
#rm $outputsPath"/"$runNum"_calls.flt-qualDP.bcf"
#rm $outputsPath"/"$runNum"_calls.flt-qualDP-homo.bcf"

# status message
echo "Analysis conplete!"

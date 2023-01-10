#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_VC_jobOutput
#$ -pe smp 8
#Script to run the SSR pipeline
#Usage: qsub variantCalling_bcftools.sh inputsFile
#Usage Ex: qsub variantCalling_bcftools.sh inputPaths_romero_run1.txt
#Usage Ex: qsub variantCalling_bcftools.sh inputPaths_romero_run2.txt
#Usage Ex: qsub variantCalling_bcftools.sh inputPaths_romero_run3.txt
#Usage Ex: qsub variantCalling_bcftools.sh inputPaths_romero_run4.txt
#Usage Ex: qsub variantCalling_bcftools.sh inputPaths_romero_run5.txt

#Required modules for ND CRC servers
module load bio

#Activate the python2 environment for local run
source /afs/crc.nd.edu/user/e/ebrooks5/.bashrc
conda activate /afs/crc.nd.edu/user/e/ebrooks5/.conda/envs/python2

#Activate the python2 environment for job run
#conda activate python2

#Make sure numpy is installed
#pip install numpy

#Retrieve input argument of a inputs file
inputsFile=$1

#Retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
#Retrieve genome reference absolute path for alignment
ref=$(grep "genomeReference:" ../"InputData/"$inputsFile | tr -d " " | sed "s/genomeReference://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

# setup the inputs path
inputsPath=$outputsPath"/"$projectDir"_SSR_prep"

#Name output file of inputs
inputOutFile=$inputsPath"/pipeline_VC_summary.txt"
#Add pipeline info to outputs
echo -e "SSR pipline variant calling inputs for $projectDir \n" > $inputOutFile


#Variant Calling Stage - SNP Calling Workflow

#Set input paths
inputsPath=$inputsPath"/aligned"

#Loop through all filtered sam files
for f1 in "$inputsPath"/*filter50.sam; do
	echo "Processing file $f"
	path=$(cat $f | sed 's/\.sam//g')
	#Calculate the read coverage of positions in the genome
	bcftools mpileup --threads 8 -d 8000 -Q 20 -Ob -o "$path"_raw.bcf -f "$ref" "$f" 
	echo bcftools mpileup --threads 8 -d 8000 -Q 20 -Ob -o "$path"_raw.bcf -f "$ref" "$f" >> "$inputOutFile"
	#Detect the single nucleotide polymorphisms 
	bcftools call --threads 8 -mv -Oz -o "$path"_calls.vcf.gz "$path"_raw.bcf 
	echo bcftools call --threads 8 -mv -Oz -o "$path"_calls.vcf.gz "$path"_raw.bcf >> "$inputOutFile"
	#Index vcf file
	bcftools index --threads 8 "$path"_calls.vcf.gz
	echo bcftools index --threads 8 "$path"_calls.vcf.gz >> "$inputOutFile"
	#Normalize indels
	bcftools norm --threads 8 -f "$ref" "$path"_calls.vcf.gz -Ob -o "$path"_calls.norm.bcf
	echo bcftools norm --threads 8 -f "$ref" "$path"_calls.vcf.gz -Ob -o "$path"_calls.norm.bcf >> "$inputOutFile"
	#Filter adjacent indels within 5bp
	bcftools filter --threads 8 --IndelGap 5 "$path"_calls.norm.bcf -Ob -o "$path"_calls.norm.flt-indels.bcf
	echo bcftools filter --threads 8 --IndelGap 5 "$path"_calls.norm.bcf -Ob -o "$path"_calls.norm.flt-indels.bcf >> "$inputOutFile"
	#Include sites where FILTER is true
	#bcftools query -i'FILTER="."' -f'%CHROM %POS %FILTER\n' "$path"_calls.norm.flt-indels.bcf > "$path"_filtered.bcf
	#echo bcftools query -i'FILTER="."' -f'%CHROM %POS %FILTER\n' "$path"_calls.norm.flt-indels.bcf ">" "$path"_filtered.bcf >> "$inputOutFile"
done

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

# setup the variant calling directory
dataPath=$inputsPath"/variants"
# create the directory
mkdir $dataPath
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

#Outputs directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_SNP"

# copy software prep summary
cp $inputsPath"/software_prep_summary.txt" $outputsPath

#Name output file of inputs
versionFile=$outputsPath"/software_VC_summary.txt"
#Add pipeline info to outputs
echo -e "SSR pipeline variant calling software versions for $projectDir \n" > $versionFile
# retrieve software version
bcftools --version >> $versionFile


#Variant Calling Stage - SNP Calling Workflow

#Set input paths
inputsPath=$inputsPath"/sorted"

#Loop through all filtered sam files
for f in "$inputsPath"/*sortedCoordinate.bam; do
	# remove the file extension
	path=$(echo $dataPath"/"$fileName | sed 's/\.bam$//g')
	# status message
	echo "Processing file "$path".bam ..."
	#Calculate the read coverage of positions in the genome
	bcftools mpileup --threads 8 -d 8000 -Ob -o $path"_raw.bcf" -f $ref $f 
	#Detect the single nucleotide polymorphisms 
	bcftools call --threads 8 -mv -Oz -o $path"_calls.vcf.gz" $path"_raw.bcf" 
	#Index vcf file
	bcftools index --threads 8 $path"_calls.vcf.gz"
	#Normalize indels
	bcftools norm --threads 8 -f $ref $path"_calls.vcf.gz" -Ob -o $path"_calls.norm.bcf"
	#Filter adjacent indels within 5bp
	bcftools filter --threads 8 --IndelGap 5 $path"_calls.norm.bcf" -Ob -o $path"_calls.norm.flt-indels.bcf"
	# status message
	echo "Processed!"
done

# status message
echo "Analysis conplete!"

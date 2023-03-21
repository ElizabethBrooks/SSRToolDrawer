#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N test_ssr_SNP_jobOutput
#$ -pe smp 4

# script to run the SSR pipeline
# usage: qsub ssr_pipeline_SNP_test.sh runInputs
# usage Ex: qsub ssr_pipeline_SNP_test.sh inputs_run1.txt
# usage Ex: qsub ssr_pipeline_SNP_test.sh inputs_run2.txt
# usage Ex: qsub ssr_pipeline_SNP_test.sh inputs_run3.txt
# usage Ex: qsub ssr_pipeline_SNP_test.sh inputs_run4.txt
# usage Ex: qsub ssr_pipeline_SNP_test.sh inputs_run5.txt

# required modules for ND CRC servers
module load bio
module load parallel

# activate the python2 environment for local run
source /afs/crc.nd.edu/user/e/ebrooks5/.bashrc
conda activate /afs/crc.nd.edu/user/e/ebrooks5/.conda/envs/python2

## if necessary
## activate the python2 environment for job run
##conda activate python2
## make sure numpy is installed
##pip install numpy

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve current working directory
currDir=$(pwd)

# retrieve base directory path
baseDir=$(dirname $currDir)

# retrieve the run number 
runNum=$(grep "run:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/run://g")
# retrieve the project ID 
projectDir=$(grep "ID:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
# retrieve ssr info path
infoPath=$(grep "ssrInfo:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrInfo://g")
# retrieve ssr regions path
regionsPath=$(grep "ssrRegions:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrRegions://g")
# retrieve primers path
primerPath=$(grep "primers:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/primers://g")
# bamclipper tool path
clipperPath=$(grep "bamclipperTool:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/bamclipperTool://g")
# retrieve reference path
ref=$(grep "reference:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")
# retrieve analysis outputs path
outputsPath=$(grep "outputs:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_SNP"

# setup the inputs path
inputsPath=$outputsPath"/"$projectDir"_SSR_prep"

# SSR Analysis Stage - SNP Calling Workflow

# status message
echo "SSR SNP analysis started..."

# move to the inputs directory
cd $currDir"/Scripts"

# perform sorting and variant calling
#bash sorting_samtools.sh $inputsPath $projectDir
#bash variantCalling_bcftools.sh $inputsPath $projectDir $ref $runNum

# remove header lines from the vcf file
for f2 in $inputsPath"/variants/"*"_calls.norm.bcf"; do
	# print status message
	echo "Removing header from $f2"
	# remove extension
	newName=$(echo $f2 | sed 's/\.vcf//g')
	# convert bcf to vcf
	bcftools convert -Ov -o $newName".vcf" $f2
	# remove header lines
	grep -v "#" $newName".vcf" > $newName".noHeader.vcf"
	# status message
	echo "Processed!"
done

# retrieve and format sample tag list
#sampleTags=$(for i in $inputsPath"/variants/"*".noHeader.vcf"; do basename $i | sed "s/^/\"/g" | sed "s/\.noHeader\.vcf$/\",/g" | tr '\n' ' '; done)
#sampleTags=$(echo $sampleTags | sed 's/.$//')

# find and replace the sample list
#sed -i "s/\"FIND_ME_REPLACE_ME\"/$sampleTags/g" Format_VCF-Matrix.py

# format matrix
#python2 Format_VCF-Matrix.py

# re-name and move output matrix
#mv VCF_Matrix.txt $outputsPath"/"$runNum".txt"

# status message
echo "SSR VC analysis complete!"

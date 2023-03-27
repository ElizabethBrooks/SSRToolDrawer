#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N test_ssr_SNP_jobOutput
#$ -pe smp 4

# script to run the SSR pipeline
# usage: qsub ssr_pipeline_SNP.sh runInputs
# usage Ex: qsub ssr_pipeline_SNP_test.sh inputs_run1.txt

# TO-DO
# add the combination of plate runs to VC

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
# retrieve reference path
ref=$(grep "reference:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")
# retrieve analysis outputs path
outputsPath=$(grep "outputs:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_SNP"
mkdir $outputsPath
# check if the folder already exists
#if [ $? -ne 0 ]; then
#	echo "The $outputsPath directory already exsists... please remove before proceeding."
#	exit 1
#fi

# setup the inputs path
inputsPath=$outputsPath"/"$projectDir"_SSR_prep"
mkdir $inputsPath


# SSR Analysis Stage - SNP Calling Workflow

# status message
echo "SSR SNP analysis started..."

# move to the alignment directory
cd $inputsPath"/aligned"

# copy pipeline scripts to the aligned directory
cp $baseDir"/SNP_Calling/Scripts/SamIAm.py" $inputsPath"/aligned"



# move to pipeline scripts directory
cd $currDir"/Scripts"


# run script to clip primer and ssr sequences
bash clipping_samtools_bamclipper.sh $inputsPath $baseDir

# move to pipeline scripts directory
cd $currDir"/Scripts"

# TO-DO
# consider filtering by mapping quality

# consider running script to perform variant calling
bash variantCalling_bcftools.sh $inputsPath $inputsFile $baseDir

# move to variants directory
cd $inputsPath"/variants"

# subset vcf file by sample and remove header lines
for f2 in $inputsPath"/clipped/"*".readGroups.bam"; do
	# retrieve sample name and remove the file extension
	sampleTag=$(basename $f2 | sed 's/\.readGroups\.bam$//g')
	# print status message
	echo "Subsetting VCF and removing header for $sampleTag"
	# subset vcf files by sample and remove header
	bcftools view --threads 4 -H -Ov -o $inputsPath"/variants/"$sampleTag".noHeader.vcf" -s $sampleTag $inputsPath"/variants/"$runNum"_noSSR.vcf"	
	# status message
	echo "Processed!"
done

# retrieve and format sample tag list
sampleTagList=$(for i in $inputsPath"/variants/"*".noHeader.vcf"; do basename $i | sed "s/^/\"/g" | sed "s/\.noHeader\.vcf$/\",/g" | tr '\n' ' '; done)
sampleTagList=$(echo $sampleTagList | sed 's/.$//')

# copy pipeline scripts to the variants directory
cp $currDir"/Scripts/Format_VCF-Matrix.py" $inputsPath"/variants"

# find and replace the sample list
sed -i "s/\"FIND_ME_REPLACE_ME\"/$sampleTagList/g" $inputsPath"/variants/"Format_VCF-Matrix.py

# format matrix
python2 Format_VCF-Matrix.py

# re-name and move output matrix
mv VCF_Matrix.txt $outputsPath"/"$runNum".txt"

# clean up
#rm -r $inputsPath

# status message
echo "SSR VC analysis complete!"



#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N test_ssr_SNP_jobOutput
#$ -pe smp 4

# script to run the SSR pipeline
# usage: qsub ssr_pipeline_SNP.sh runList
# usage Ex: qsub ssr_pipeline_SNP_test.sh run1 run2 run3
# usage Ex: qsub ssr_pipeline_SNP_test.sh run4 run5 run6 run7
# usage Ex: qsub ssr_pipeline_SNP_test.sh run1 run2 run3 run4 run5 run6 run7

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve current working directory
currDir=$(pwd)

# retrieve base directory path
baseDir=$(dirname $currDir)

# retrieve analysis outputs path
outputsPath=$(grep "outputs:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for project analysis
outputsPath=$outputsPath"/SSR_SNP"
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "Warning! The $outputsPath directory already exsists. Proceeding..."
fi

# SSR Analysis Stage - SNP Calling Workflow

# status message
echo "SSR SNP analysis started..."

# move to scripts dorectory
cd $baseDir"/SNP_Calling/Scripts"

# loop over each input file
for run in "$@"; do
	# status message
	echo "Preparing data for $run"
	# set inputs file name
	subsetFile="inputs_"$run".txt"
	# run bash script to process the current subset
	#qsub ssr_pipeline_subset_SNP.sh $subsetFile $baseDir
	bash ssr_pipeline_subset_SNP.sh $subsetFile $baseDir
	# setup the inputs path
	#subsetPath=$outputsPath"/SSR_SNP_prep_"$run
done

# load software modules
#module load bio/2.0

# setup the inputs path
#inputsPath=$outputsPath"/SSR_SNP_prep"

## TO-DO: combine subset prep files

# run script to perform variant calling
#bash variantCalling_bcftools.sh $inputsPath $baseDir $inputsPathList

# run script to perform variant filtering
#bash variantFiltering_bcftools.sh $inputsPath $baseDir

# run script to perform variant trimming
#bash variantTrimming_bedtools.sh $inputsPath $baseDir

# format matrix
#bash variantMatrix_bcftools.sh $inputsPath $baseDir

# clean up
#rm -r $inputsPath

# status message
echo "SSR VC analysis complete!"

#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_SNP_jobOutput
#$ -q largemem
#$ -pe smp 8

# script to run the SSR pipeline
# usage: qsub ssr_pipeline_SNP.sh runList
# usage Ex: qsub ssr_pipeline_SNP.sh run1 run2 run3 run4 run5 run6 run7 run8
## job 795810
# usage Ex: qsub ssr_pipeline_SNP.sh run1 run2 run3 run4 run5 run6 run7 run8 run9
## job 795812
# usage Ex: qsub ssr_pipeline_SNP.sh run1 run2 run3 run4 run5 run6 run7 run8 run9 run10
## job 783901 -> SNP_Calling_run1_to_run10_test2
## job 792141 -> SNP_Calling_run1_to_run10_test3
## job 793827 -> SNP_Calling_run1_to_run10_test4
# usage Ex: qsub ssr_pipeline_SNP.sh run1 run2 run3 run4 run5 run6 run7 run8 run10
## job 

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve current working directory
currDir=$(pwd)

# retrieve base directory path
baseDir=$(dirname $currDir)

# retrieve analysis outputs path
outputsPath=$(grep "outputs:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for project analysis
#outputsPath=$outputsPath"/SNP_Calling"
outputsPath=$outputsPath"/SNP_Calling_"$1"_to_"${@: -1}
mkdir $outputsPath

# comment this out to re-start the analysis
# before re-starting the analysis, make sure to remove any sub directories that were not completely analyzed
# this should be the last directory that was created during the analysis
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# make a directory for the software version info
mkdir $outputsPath"/info"

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
	bash ssr_pipeline_subset_SNP.sh $subsetFile $baseDir $outputsPath
done

# load software modules
module load bio/2.0

# run script to perform variant calling
bash variantCalling_bcftools.sh $baseDir $outputsPath

# run script to perform variant filtering
bash variantFiltering_bcftools.sh $baseDir $outputsPath

# run script to perform variant trimming
bash variantTrimming_bedtools.sh $baseDir $outputsPath

# format matrix
bash variantMatrix_bcftools.sh $baseDir $outputsPath

# clean up
#rm -r $outputsPath"/SNP_Calling_prep_run"*
#rm -r $outputsPath"/variantsCalled"
#rm -r $outputsPath"/variantsFiltered"
#rm -r $outputsPath"/variantsTrimmed"

# output run time
echo $SECONDS > $outputsPath"/combined_seconds.txt"

# status message
echo "SSR VC analysis complete!"

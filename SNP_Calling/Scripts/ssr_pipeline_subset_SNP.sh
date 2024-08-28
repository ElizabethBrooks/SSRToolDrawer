#!/bin/bash

# script to run the SSR pipeline
# usage: bash ssr_pipeline_SNP.sh inputsFile baseDir

# load required software modules
module load bio/2.0
module load parallel

## if necessary
## activate the python2 environment for job run
##conda activate python2
## make sure numpy is installed
##pip install numpy

# retrieve input argument of a inputs file
inputsFile=$1

# retrieve base directory path
baseDir=$2

# retrieve the run number 
runNum=$(grep "run:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/run://g")
# retrieve ssr info path
infoPath=$(grep "ssrInfo:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrInfo://g")
# retrieve reference path
ref=$(grep "reference:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")
# retrieve analysis outputs path
outputsPath=$(grep "outputs:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for project analysis
outputsPath=$outputsPath"/SNP_Calling"

# setup the inputs path
inputsPath=$outputsPath"/SNP_Calling_prep_"$runNum
mkdir $inputsPath

# before re-starting the analysis, make sure to remove any sub directories that were not completely analyzed
# this should be the last directory that was created during the analysis
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "Warning! The $inputsPath directory already exsists. Proceeding..."
fi

# prepare data for analysis
cd $baseDir"/Prep"
bash ssr_pipeline_prep.sh $inputsPath $inputsFile $baseDir

## TO-DO
## make sure to check mapping efficiency
## note that we did not remove primers in advance
## consider demultiplexing using primers
## or removing primers using trimmomatic along with adapters

## TO-DO
## consider sorting and removinng PCR duplicates before filtering with SamIAM.py


# SSR Analysis Stage - SNP Calling Workflow

# unload modules
module unload bio/2.0
module unload parallel

# activate the python2 environment for local run
source /afs/crc.nd.edu/user/e/ebrooks5/.bashrc
conda activate /afs/crc.nd.edu/user/e/ebrooks5/.conda/envs/python2

# status message
echo "SSR SNP data prep started for $runNum..."

# move to the alignment directory
cd $inputsPath"/aligned"

# copy pipeline scripts to the aligned directory
cp $baseDir"/SNP_Calling/Scripts/SamIAm.py" $inputsPath"/aligned"

# set outputs path
outputsDir=$inputsPath"/alignedFiltered50"
# create the directory
mkdir $outputsDir

# loop through all aligned sam files
for f1 in $inputsPath"/aligned/"*".sam"; do
	# trim file path from current folder name
	curSampleNoPath=$(basename "$f1" | sed 's/\.sam$//g')
	# print status message
	echo "Processing $f1"
	# run SSR pipeline
	python2 SamIAm.py -sam $f1 -C $infoPath -p "yes"
	# replace SAM header
	grep "^@" $f1 > $outputsDir"/"$curSampleNoPath".header.sam"
	# append filtered sequences
	cat $inputsPath"/aligned/"$curSampleNoPath".sam.filter50.sam" >> $outputsDir"/"$curSampleNoPath".header.sam"
	rm $inputsPath"/aligned/"$curSampleNoPath".sam.filter50.sam"
	rm $inputsPath"/aligned/"$curSampleNoPath".sam.hitInfo"
done

# deactivate the conda environment
conda deactivate

# load required software modules
module load bio/2.0
module load parallel

# move to pipeline scripts directory
cd $baseDir"/SNP_Calling/Scripts"

# run script to perform sorting and removal of pcr duplicates
bash sorting_samtools.sh $inputsPath $baseDir

# run script to keep only unique read alignments
bash filterByMapQ_samtools.sh $inputsPath $baseDir

# run script to clip primer and ssr sequences
bash clipping_samtools_bamclipper.sh $inputsPath $baseDir $runNum

## TO-DO
## consider single sample calling
## https://github.com/samtools/bcftools/issues/811
## bcftools mpileup -a AD
## bcftools call -G -

# clean up
#rm -r $inputsPath

# status message
echo "SSR VC data prep complete for $runNum!"

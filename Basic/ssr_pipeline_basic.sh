#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_basic_jobOutput
#$ -pe smp 8

# script to run the SSR pipeline
# usage: qsub ssr_pipeline_basic.sh inputsFile
# usage Ex: qsub ssr_pipeline_basic.sh inputPaths_romero_run1.txt
# usage Ex: qsub ssr_pipeline_basic.sh inputPaths_romero_run2.txt
# usage Ex: qsub ssr_pipeline_basic.sh inputPaths_romero_run3.txt
# usage Ex: qsub ssr_pipeline_basic.sh inputPaths_romero_run4.txt
# usage Ex: qsub ssr_pipeline_basic.sh inputPaths_romero_run5.txt
# usage Ex: qsub ssr_pipeline_basic.sh inputPaths_romero_run6.txt

# required modules for ND CRC servers
module load bio

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

# retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
# retrieve adapter absolute path for alignment
infoPath=$(grep "info:" ../"InputData/"$inputsFile | tr -d " " | sed "s/info://g")
# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

# make a new directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_basic"
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# setup the inputs path
inputsPath=$outputsPath"/"$projectDir"_SSR_prep"
mkdir $inputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $inputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# prepare data for analysis
cd ../Prep
bash ssr_pipeline_prep.sh $inputsFile $inputsPath


# SSR Analysis Stage - Basic Workflow

# status message
echo "SSR basic analysis started..."

# copy pipeline scripts to inputs directory
cd ../Basic
cp GapGenes.v3.py $inputsPath"/aligned"
cp SnipMatrix.py $inputsPath"/aligned"
cp Format_Matrix.py $inputsPath"/aligned"

# move to the inputs directory
cd $inputsPath"/aligned"

# loop through all filtered sam files
for f1 in $inputsPath"/aligned/"*".sam"; do
	# print status message
	echo "Processing $f1"
	# run SSR pipeline
	python2 GapGenes.v3.py -sam $f1 -C $infoPath -P "unpaired"
	python2 SnipMatrix.py $f1".Matrix.txt"
	# status message
	echo "Processed!"
done

# retrieve and format sample tag list
sampleTags=$(for i in $inputsPath"/aligned/"*.sam; do basename $i | sed "s/^/\"/g" | sed "s/\.sam/\",/g" | tr '\n' ' '; done)
sampleTags=$(echo $sampleTags | sed 's/.$//')

# find and replace the sample list
sed -i "s/\"FIND_ME_REPLACE_ME\"/$sampleTags/g" Format_Matrix.py

# format matrix
python2 Format_Matrix.py

# re-name and move output matrix
mv SNP_Matrix.txt $outputsPath"/"$projectDir"_SNP_Matrix.txt"

# clean up
rm -r $inputsPath

# status message
echo "SSR SNP analysis complete!"

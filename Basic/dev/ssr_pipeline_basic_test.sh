#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_basic_test_jobOutput
#$ -pe smp 4

# script to run the SSR pipeline
# usage: qsub ssr_pipeline_basic_test.sh runInputs
# usage Ex: qsub ssr_pipeline_basic_test.sh inputs_run7.txt

# required modules for ND CRC servers
#module load bio/2.0

# activate the python2 environment for local run
#source /afs/crc.nd.edu/user/e/ebrooks5/.bashrc
#conda activate /afs/crc.nd.edu/user/e/ebrooks5/.conda/envs/python2

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
# retrieve analysis outputs path
outputsPath=$(grep "outputs:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_Basic"
#mkdir $outputsPath
# check if the folder already exists
#if [ $? -ne 0 ]; then
#	echo "The $outputsPath directory already exsists... please remove before proceeding."
#	exit 1
#fi

# setup the inputs path
inputsPath=$outputsPath"/"$projectDir"_SSR_prep"
#mkdir $inputsPath

# prepare data for analysis
#cd $baseDir"/Prep"
#bash ssr_pipeline_prep.sh $inputsPath $inputsFile $baseDir


# SSR Analysis Stage - Basic Workflow

# unload modules
#module unload bio/2.0

# activate the python2 environment for local run
source /afs/crc.nd.edu/user/e/ebrooks5/.bashrc
conda activate /afs/crc.nd.edu/user/e/ebrooks5/.conda/envs/python2

# status message
echo "SSR basic analysis started..."

# copy pipeline scripts to inputs directory
cp $baseDir"/Basic/Scripts/"* $inputsPath"/aligned"

# move to the inputs directory
cd $inputsPath"/aligned"

# loop through all filtered sam files
#for f1 in $inputsPath"/aligned/"*".sam"; do
	# print status message
	#echo "Processing $f1"
	# run SSR pipeline
	#python2 GapGenes.v3.py -sam $f1 -C $infoPath
	##python2 GapGenes.v3.py -sam $f1 -C $infoPath -P "paired"
	#python2 SnipMatrix.py $f1".Matrix.txt"
	##python2 SnipMatrix.py $f1".Matrix.1.txt"
	##python2 SnipMatrix.py $f1".Matrix.2.txt"
	# status message
	#echo "Processed!"
#done

# retrieve and format sample tag list
sampleTags=$(for i in $inputsPath"/aligned/"*.sam; do basename $i | sed "s/^/\"/g" | sed "s/\.sam/\",/g" | tr '\n' ' '; done)
sampleTags=$(echo $sampleTags | sed 's/.$//')

# find and replace the sample list
sed -i "s/\"FIND_ME_REPLACE_ME\"/$sampleTags/g" Format_Matrix.py

# format matrix
python2 Format_Matrix.py

# re-name and move output matrix
mv SNP_Matrix.txt $outputsPath"/"$runNum".txt"

# clean up
#rm -r $inputsPath

# status message
echo "SSR SNP analysis complete!"

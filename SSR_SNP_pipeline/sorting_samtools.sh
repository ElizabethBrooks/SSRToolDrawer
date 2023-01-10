#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_VC_jobOutput
#$ -pe smp 8
#Script to run the SSR pipeline
#Usage: qsub sorting_samtools.sh inputsFile
#Usage Ex: qsub sorting_samtools.sh inputPaths_romero_run1.txt
#Usage Ex: qsub sorting_samtools.sh inputPaths_romero_run2.txt
#Usage Ex: qsub sorting_samtools.sh inputPaths_romero_run3.txt
#Usage Ex: qsub sorting_samtools.sh inputPaths_romero_run4.txt
#Usage Ex: qsub sorting_samtools.sh inputPaths_romero_run5.txt

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
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

# setup the inputs path
inputsPath=$outputsPath"/"$projectDir"_SSR_prep"

# setup the variant calling directory
dataPath=$inputsPath"/sorted"
# create the directory
mkdir $dataPath
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

#Outputs directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_SNP"

#Name output file of inputs
versionFile=$outputsPath"/software_sorting_summary.txt"
#Add pipeline info to outputs
echo -e "SSR pipeline sorting software versions for $projectDir \n" > $versionFile
# retrieve software version
samtools --version >> $versionFile


#Sorting Stage - SNP Calling Workflow

#Set input paths
inputsPath=$inputsPath"/aligned"

#Loop through all filtered sam files
for f in "$inputsPath"/*filter50.sam; do
	# remove two file extensions
	pathSam=$(echo $f | sed 's/\.filter50\.sam$//g')
	# remove the path from the file name
	fileName=$(basename $f)
	# add the sam header to the input file
	grep "^@" $pathSam > $dataPath"/"$fileName
	cat $f >> $dataPath"/"$fileName
	# remove the file extension
	path=$(echo $dataPath"/"$fileName | sed 's/\.sam$//g')
	# trim file path from current folder name
	curSampleNoPath=$(basename "$f" | sed 's/\.sam$//g')
	# status message
	echo "Processing file "$path".sam ..."
	# convert output sam files to bam format for downstream analysis
	samtools view -@ 8 -bS $dataPath"/"$fileName > $path".bam"
	#Run samtools to prepare mapped reads for sorting
	#using 8 threads
	samtools sort -@ 8 -n -o $path".sortedName.bam" -T "/tmp/"$curSampleNoPath".sortedName.bam" $path".bam"
	#Run fixmate -m to update paired-end flags for singletons
	samtools fixmate -m $path".sortedName.bam" $path".sortedFixed.bam"
	#Clean up
	rm $path".sortedName.bam"
	#Run samtools to prepare mapped reads for sorting by coordinate
	#using 8 threads
	samtools sort -@ 8 -o $path".sortedCoordinate.bam" -T "/tmp/"$curSampleNoPath".sortedCoordinate.bam" $path".sortedFixed.bam"
	# clean up
	rm $path".sortedFixed.bam"
	#Remove duplicate reads
	samtools markdup -r $path".sortedCoordinate.bam" $path".noDups.bam"
	# status message
	echo "Processed!"
done

# status message
echo "Analysis conplete!"

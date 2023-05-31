#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_SNP_jobOutput
#$ -pe smp 4

# script to run the SSR pipeline
# usage: qsub ssr_pipeline_SNP.sh runInputs
# usage Ex: qsub ssr_pipeline_SNP.sh inputs_run1.txt inputs_run2.txt inputs_run3.txt inputs_run4.txt inputs_run5.txt inputs_run6.txt

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

# retrieve current working directory
currDir=$(pwd)

# retrieve base directory path
baseDir=$(dirname $currDir)

# loop over each input file
for arg in "$@"; do
	# run bash cript to process the current subset
	bash ssr_pipeline_subset_SNP.sh $arg
done

# TO-DO
# consider single sample calling
# https://github.com/samtools/bcftools/issues/811
# bcftools mpileup -a AD
# bcftools call -G -

# run script to perform variant calling
bash variantCalling_bcftools.sh $inputsPath $inputsFile $baseDir

# run script to perform variant filtering
bash variantFiltering_bcftools.sh $inputsPath $inputsFile $baseDir

# run script to perform variant trimming
bash variantTrimming_bedtools.sh $inputsPath $inputsFile $baseDir

# format matrix
bash variantMatrix_bcftools.sh $inputsPath $inputsFile $baseDir

# clean up
#rm -r $inputsPath

# status message
echo "SSR VC analysis complete!"

#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_SNP_jobOutput
#$ -pe smp 4

# script to run the SSR pipeline
# usage: qsub ssr_pipeline_SNP_test.sh inputsFile
# usage Ex: qsub ssr_pipeline_SNP_test.sh inputPaths_romero_run1.txt
# usage Ex: qsub ssr_pipeline_SNP_test.sh inputPaths_romero_run2.txt
# usage Ex: qsub ssr_pipeline_SNP_test.sh inputPaths_romero_run3.txt
# usage Ex: qsub ssr_pipeline_SNP_test.sh inputPaths_romero_run4.txt
# usage Ex: qsub ssr_pipeline_SNP_test.sh inputPaths_romero_run5.txt

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
outputsPath=$outputsPath"/"$projectDir"_SSR_SNP_test"
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
bash ssr_pipeline_prep.sh $inputsFile $outputsPath


# SSR Analysis Stage - SNP Calling Workflow

# status message
echo "SSR SNP analysis started..."

# copy pipeline scripts to inputs directory
cd ../SNP_Calling
cp SamIAm.py $inputsPath"/aligned"
cp sorting_samtools.sh $inputsPath"/aligned"
cp variantCalling_bcftools.sh $inputsPath"/aligned"
cp Format_VCF-Matrix.py $inputsPath"/variants"

# move to the inputs directory
cd $inputsPath"/aligned"

# loop through all filtered sam files
for f1 in $inputsPath"/aligned/"*".sam"; do
	# print status message
	echo "Processing $f1"
	# run SSR pipeline
	python2 SamIAm.py -sam $f1 -C $infoPath -p "yes"
	# status message
	echo "Processed!"
done

# perform sorting and variant calling
bash sorting_samtools.sh $inputsFile $inputsPath
bash variantCalling_bcftools.sh $inputsFile $inputsPath

# move to the inputs directory
cd $inputsPath"/variants"

# remove headers from the vcf files
for f2 in $inputsPath"/variants/"*".flt-indels.vcf"; do
	# print status message
	echo "Removing header from $f2"
	# create new file name
	newName=$(echo $f2 | sed 's/\.sam\.filter50\.sortedCoordinate\_calls\.norm\.flt\-indels\.vcf/\.noHeader\.vcf/g')
	# remove header
	grep -v "#" $f2 > $newName
	# status message
	echo "Processed!"
done

# retrieve and format sample tag list
sampleTags=$(for i in $inputsPath"/variants/"*".noHeader.vcf"; do basename $i | sed "s/^/\"/g" | sed "s/\.noHeader\.vcf$/\",/g" | tr '\n' ' '; done)
sampleTags=$(echo $sampleTags | sed 's/.$//')

# find and replace the sample list
sed -i "s/\"FIND_ME_REPLACE_ME\"/$sampleTags/g" Format_VCF-Matrix.py

# format matrix
python2 Format_VCF-Matrix.py

# re-name and move output matrix
mv VCF_Matrix.txt $outputsPath"/"$projectDir"_VCF_Matrix.txt"

# clean up
#rm -r $inputsPath

# status message
echo "SSR VC analysis complete!"

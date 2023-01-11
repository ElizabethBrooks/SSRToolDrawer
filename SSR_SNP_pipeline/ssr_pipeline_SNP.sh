#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_SNP_jobOutput
#Script to run the SSR pipeline
#Usage: qsub ssr_pipeline_SNP.sh inputsFile
#Usage Ex: qsub ssr_pipeline_SNP.sh inputPaths_romero_run1.txt
#Usage Ex: qsub ssr_pipeline_SNP.sh inputPaths_romero_run2.txt
#Usage Ex: qsub ssr_pipeline_SNP.sh inputPaths_romero_run3.txt
#Usage Ex: qsub ssr_pipeline_SNP.sh inputPaths_romero_run4.txt
#Usage Ex: qsub ssr_pipeline_SNP.sh inputPaths_romero_run5.txt

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
#Retrieve adapter absolute path for alignment
infoPath=$(grep "info:" ../"InputData/"$inputsFile | tr -d " " | sed "s/info://g")
#Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"InputData/"$inputsFile | tr -d " " | sed "s/outputs://g")

# setup the inputs path
inputsPath=$outputsPath"/"$projectDir"_SSR_prep"

#Make a new directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_SNP"
mkdir $outputsPath
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi


#SSR Analysis Stage - SNP Calling Workflow

#Copy pipeline scripts to inputs directory
#cp SamIAm.py $inputsPath"/aligned"
#cp sorting_samtools.sh $inputsPath"/aligned"
#cp variantCalling_bcftools.sh $inputsPath"/aligned"
cp Format_VCF-Matrix.py $inputsPath"/variants"

#Move to the inputs directory
#cd $inputsPath"/aligned"

#Loop through all filtered sam files
#for f in $inputsPath"/aligned/"*".sam"; do
	#Print status message
#	echo "Processing $f"
	#Run SSR pipeline
#	python2 SamIAm.py -sam $f -C $infoPath -p "yes"
	#Status message
#	echo "Processed!"
#done

# perform variant calling
#bash sorting_samtools.sh $inputsFile
#bash variantCalling_bcftools.sh $inputsFile

#Move to the inputs directory
cd $inputsPath"/variants"

#Retrieve and format sample tag list
sampleTags=$(for i in $inputsPath"/variants/"*".sam.filter50.sortedCoordinate_calls.norm.flt-indels.vcf"; do basename $i | sed "s/^/\"/g" | sed "s/\.vcf$/\",/g" | tr '\n' ' '; done)
sampleTags=$(echo $sampleTags | sed 's/.$//')

#Find and replace the sample list
sed -i "s/\"FIND_ME_REPLACE_ME\"/$sampleTags/g" Format_VCF-Matrix.py

#Format matrix
python2 Format_VCF-Matrix.py

#Clean up
#rm $inputsPath"/aligned/SamIAm.py"
#rm $inputsPath"/aligned/sorting_samtools.sh"
#rm $inputsPath"/aligned/variantCalling_bcftools.sh"
rm $inputsPath"/variants/Format_VCF-Matrix.py"
#rm -r $inputsPath

#Re-name and move output matrix
mv VCF_Matrix.txt $outputsPath"/"$projectDir"_VCF_Matrix.txt"

#Print status message
echo "Analysis complete!"

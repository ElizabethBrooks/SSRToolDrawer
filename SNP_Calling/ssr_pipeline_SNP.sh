#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N ssr_SNP_jobOutput
#$ -pe smp 4

# script to run the SSR pipeline
# usage: qsub ssr_pipeline_SNP.sh runInputs
# usage Ex: qsub ssr_pipeline_SNP.sh inputs_run1.txt
# usage Ex: qsub ssr_pipeline_SNP.sh inputs_run2.txt
# usage Ex: qsub ssr_pipeline_SNP.sh inputs_run3.txt
# usage Ex: qsub ssr_pipeline_SNP.sh inputs_run4.txt
# usage Ex: qsub ssr_pipeline_SNP.sh inputs_run5.txt

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
# retrieve ssr regions path
regionsPath=$(grep "ssrRegions:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrRegions://g")
# retrieve primers path
primerPath=$(grep "primers:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/primers://g")
# bamclipper tool path
clipperPath=$(grep "bamclipperTool:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/bamclipperTool://g")
# retrieve reference path
ref=$(grep "reference:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")
# retrieve analysis outputs path
outputsPath=$(grep "outputs:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for project analysis
outputsPath=$outputsPath"/"$projectDir"_SSR_SNP"
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# setup the inputs path
inputsPath=$outputsPath"/"$projectDir"_SSR_prep"
mkdir $inputsPath

# prepare data for analysis
cd $baseDir"/Prep"
bash ssr_pipeline_prep.sh $inputsFile $inputsPath

# TO-DO
# make sure to check mapping efficiency
# note that we did not remove primers in advance
# consider demultiplexing using primers
# or removing primers using trimmomatic along with adapters


# SSR Analysis Stage - SNP Calling Workflow

# status message
echo "SSR SNP analysis started..."

# copy pipeline scripts to alignment directory
cp $baseDir"/SNP_Calling/Scripts/SamIAm.py" $inputsPath"/aligned"
cp -r $clipperPath"/"* $inputsPath"/aligned"

# move to the alignment directory
cd $inputsPath"/aligned"

# loop through all aligned sam files
for f1 in $inputsPath"/aligned/"*".sam"; do
	# remove file extension
	noExt=$(echo $f1 | sed 's/\.sam//g')
	# remove file path
	sample=$(basename $f1 | sed 's/\.sam//g')
	# print status message
	echo "Processing $f1"
	# run SSR pipeline
	python2 SamIAm.py -sam $f1 -C $infoPath -p "yes"
	# replace SAM header
	grep "^@" $f1 > $noExt".header.sam"
	# append filtered sequences
	cat $noExt".filter50.sam" >> $noExt".header.sam"
	rm $noExt".filter50.sam"
	# convert sam to bam
	samtools view -@ 4 -bo $noExt".header.bam" $noExt".header.sam"
	rm $noExt".header.sam"
	# index the bam file
	samtools index $noExt".header.bam" 
	# remove primers sequences
	./bamclipper.sh -b $noExt".header.bam" -p $primerPath -n 4
	rm $noExt".header.bam"
	# remove SSR sequences
	samtools view -@ 4 -bo $noExt".overlap.bam" -U $noExt".noOverlap.bam" -L $regionsPath $noExt".header.primerclipped.bam"
	rm $noExt".header.primerclipped.bam"
	rm $noExt".overlap.bam"
	# add read groups
	samtools addreplacerg -@ 4 -r ID:"SSR_"$runNum"_"$sample -r SM:$sample -o $noExt".readGroups.bam" $noExt".noOverlap.bam"
	rm $noExt".noOverlap.bam"
	# status message
	echo "Processed!"
done

# TO-DO
# consider merging BAM files before variant calling

# move to pipeline scripts directory
cd $currDir"/Scripts"

# perform sorting and variant calling
bash sorting_samtools.sh $inputsPath $projectDir
bash variantCalling_bcftools.sh $inputsPath $projectDir $ref $runNum

# remove header lines from the vcf file
for f2 in $inputsPath"/variants/"*"_calls.norm.bcf"; do
	# print status message
	echo "Removing header from $f2"
	# create new file name
	newName=$(echo $f2 | sed 's/\.vcf/\.noHeader\.vcf/g')
	# remove header lines
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
mv VCF_Matrix.txt $outputsPath"/"$runNum".txt"

# clean up
rm -r $inputsPath

# status message
echo "SSR VC analysis complete!"

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

# retrieve the run number 
runNum=$(grep "run:" ../"InputData/"$inputsFile | tr -d " " | sed "s/run://g")
# retrieve the project ID 
projectDir=$(grep "ID:" ../"InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
# retrieve ssr info path
infoPath=$(grep "ssrInfo:" ../"InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrInfo://g")
# retrieve ssr regions path
regionsPath=$(grep "ssrRegions:" ../"InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrRegions://g")
# retrieve primers path
primerPath=$(grep "primers:" ../"InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/primers://g")
# retrieve analysis outputs path
outputsPath=$(grep "outputs:" ../"InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/outputs://g")
# bamclipper tool path
clipperPath=$(grep "bamclipperTool:" ../"InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/bamclipperTool://g")

# make a new directory for project analysis
inputsPath=$outputsPath"/"$projectDir"_SSR_SNP"
mkdir $inputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $inputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# setup the downstream inputs path
inputsPath=$inputsPath"/"$projectDir"_SSR_prep"
mkdir $inputsPath

# prepare data for analysis
cd ../Prep
bash ssr_pipeline_prep.sh $inputsFile $inputsPath

# TO-DO
# make sure to check mapping efficiency
# note that we did not remove primers in advance
# consider removing using trimmomatic along with adapter seqs


# SSR Analysis Stage - SNP Calling Workflow

# status message
echo "SSR SNP analysis started..."

# move to pipeline scripts directory
cd ../SNP_Calling/Scripts

# loop through all aligned sam files
for f1 in $inputsPath"/aligned/"*".sam"; do
	# remove file extension
	noExt=$(echo $f1 | sed 's/\.sam//g')
	# remove file path
	sample=$(basename $f1 | sed 's/\.sam//g')
	# print status message
	echo "Processing $f1"
	# run SSR pipeline
	#python2 SamIAm.py -sam $f1 -C $infoPath -p "yes"
	# replace SAM header
	#grep "^@" $f1 > $noExt".header.sam"
	# append filtered sequences
	#cat $noExt".filter50.sam" >> $noExt".header.sam"
	#rm $noExt".filter50.sam"
	# convert sam to bam
	#samtools view -@ 4 -bo $noExt".header.bam" $noExt".header.sam"
	#rm $noExt".header.sam"
	# remove primers sequences
	bash $clipperPath"/bamclipper.sh" -b $noExt".header.bam" -p $primerPath -n 4
	#rm $noExt".header.bam"
	# remove SSR sequences
	#samtools view -@ 4 -bo $noExt".overlap.bam" -U $noExt".noOverlap.bam" -L $regionsPath $noExt".header.primerclipped.bam"
	#rm $noExt".header.primerclipped.bam"
	#rm $noExt".overlap.bam"
	# add read groups
	#samtools addreplacerg -@ 4 -r ID:SSR_$runNum_$sample -r SM:$sample -o $noExt".RG.bam" $noExt".noOverlap.bam"
	#rm $noExt".noOverlap.bam"
	# status message
	echo "Processed!"
done

# TO-DO
# merge BAM files before calling

# perform sorting and variant calling
#bash sorting_samtools.sh $inputsFile $outputsPath
#bash variantCalling_bcftools.sh $inputsFile $outputsPath

# TO-DO
# modify/remove this section after merging
# remove headers from the vcf files
#for f2 in $outputsPath"/variants/"*".flt-indels.vcf"; do
	# print status message
#	echo "Removing header from $f2"
	# create new file name
#	newName=$(echo $f2 | sed 's/\.sam\.filter50\.sortedCoordinate\_calls\.norm\.flt\-indels\.vcf/\.RG\.vcf/g')
	# remove header
#	grep -v "#" $f2 > $newName
	# status message
#	echo "Processed!"
#done

# retrieve and format sample tag list
#sampleTags=$(for i in $outputsPath"/variants/"*".RG.vcf"; do basename $i | sed "s/^/\"/g" | sed "s/\.RG\.vcf$/\",/g" | tr '\n' ' '; done)
#sampleTags=$(echo $sampleTags | sed 's/.$//')

# find and replace the sample list
#sed -i "s/\"FIND_ME_REPLACE_ME\"/$sampleTags/g" Format_VCF-Matrix.py

# format matrix
#python2 Format_VCF-Matrix.py

# re-name and move output matrix
#mv VCF_Matrix.txt $outputsPath"/"$runNum".txt"

# clean up
#rm -r $outputsPath"/"$projectDir"_SSR_prep"

# status message
echo "SSR VC analysis complete!"

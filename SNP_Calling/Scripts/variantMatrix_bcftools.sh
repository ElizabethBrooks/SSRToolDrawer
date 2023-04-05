#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N SNP_matrix_jobOutput

# script to create a formatted variant matrix
# usage: bash variantMatrix_bcftools.sh inputsPath inputsFile baseDir

# retrieve input outputs path
inputsPath=$1

# retrieve inputs file
inputsFile=$2

# retrieve base of working directory
baseDir=$3

# retrieve the run number 
runNum=$(grep "run:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/run://g")
# retrieve the project ID 
projectDir=$(grep "ID:" $baseDir"/InputData/"$inputsFile | tr -d " " | sed "s/ID://g")
# retrieve analysis outputs path
outputsPath=$(grep "outputs:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/outputs://g")
# retrieve ssr regions path
regionsPath=$(grep "ssrRegions:" $baseDir"/InputData/inputs_ssr_pipeline.txt" | tr -d " " | sed "s/ssrRegions://g")

# set outputs path
outputsPath=$outputsPath"/"$projectDir"_SSR_SNP"

# name of output file of inputs
versionFile=$inputsPath"/software_VC_summary.txt"

# output software version
echo "Variant matrix formatting: " >> $versionFile
bcftools --version >> $versionFile


# Variant Matrix Formatting Stage - SNP Calling Workflow

# status message
echo "Performing variant matrix formatting for $runNum"

# subset vcf file by sample and remove header lines
#for f1 in $inputsPath"/clipped/"*".readGroups.bam"; do
	# retrieve sample name and remove the file extension
#	sampleTag=$(basename $f1 | sed 's/\.readGroups\.bam$//g')
	# print status message
#	echo "Subsetting VCF and removing header for $sampleTag"
	# subset vcf files by sample and remove header
#	bcftools view --threads 4 -H -Ov -o $inputsPath"/variantsTrimmed/"$sampleTag".noHeader.vcf" -s $sampleTag $inputsPath"/variantsTrimmed/"$runNum"_trimmed.vcf"	
	# status message
#	echo "Processed!"
#done

# set output matrix file
resultsFile=$outputsPath"/"$runNum".txt"

# retrieve contig (marker) ID list
contigList=$(cat $regionsPath | cut -f 1 | tr ' ' '\t')

# add header to matrix results file
echo -en 'Sample\t'$contigList > $resultsFile

# loop over each sample
for f2 in $inputsPath"/variantsTrimmed/"*".noHeader.vcf"; do
	# retrieve sample tag
	sampleTag=$(basename $f2 | sed "s/\.noHeader\.vcf$//g")
	# output status message
	echo "Processing $sampleTag"
	# add sample tag to matrix row
	echo -en '\n'$sampleTag >> $resultsFile
	# loop over each marker
	for i in $contigList; do
		# create file with current marker variants
		cat $f2 | grep $i > $f2"."$i".tmp.txt"
		# initialize alleles
		firstAlleles="NULL"
		secondAlleles="NULL"
		# initialize line counter
		lineCount=0
		# loop over each line in the sample vcf
		while read line; do
			# increment line counter
			lineCount=$(($lineCount+1))
			# retrieve the CHROM
			chrom=$(echo $line | cut -d " " -f 1)
			# retrieve the POS
			pos=$(echo $line | cut -d " " -f 2)
			# retrieve the REF
			ref=$(echo $line | cut -d " " -f 4)
			# retrieve the ALT
			alt=$(echo $line | cut -d " " -f 5)
			# retrieve the GT alleles and translate encodings
			firstGT=$(echo $line | cut -d " " -f 10 | cut -d ":" -f 1 | cut -d "/" -f 1 | sed "s/0/$ref/g" | sed "s/1/$alt/g" | sed "s/\./NULL/g")
			secondGT=$(echo $line | cut -d " " -f 10 | cut -d ":" -f 1 | cut -d "/" -f 2 | sed "s/0/$ref/g" | sed "s/1/$alt/g" | sed "s/\./NULL/g")
			# add POS to GT alleles
			firstGT=$firstGT"("$pos")"
			secondGT=$secondGT"("$pos")"
			# check the line count
			if [[ $lineCount -eq 1 ]]; then
				# begin new list of allele variants
				firstAlleles=$firstGT
				secondAlleles=$secondGT
			else
				# append GT to lists of allele variants
				firstAlleles=$firstAlleles","$firstGT
				secondAlleles=$secondAlleles","$secondGT
			fi
		done < $f2"."$i".tmp.txt"
		# clean up
		rm $f2"."$i".tmp.txt"
		# output allele lists to the results matrix
		echo -en '\t'$firstAlleles'\t'$secondAlleles >> $resultsFile
	done < $regionsPath
done
# add a new line to the end of the matrix
echo "" >> $resultsFile

# status message
echo "Analysis conplete!"

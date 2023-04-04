#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N SNP_matrix_jobOutput

# script to create a formatted variant matrix
# usage: bash variantMatrix_formatting.sh inputsPath inputsFile baseDir


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

# set outputs path
outputsPath=$outputsPath"/"$projectDir"_SSR_SNP"

# set the inputs path
inputsPath=$inputsPath"/variantsTrimmed"

# set output matrix file
resultsFile=$outputsPath"/"$runNum".txt"


# Variant Matrix Formatting Stage - SNP Calling Workflow

# status message
echo "Performing variant matrix formatting for $runNum"

# contig (marker) ID list
contingList=$(cat $regionsPath | cut -f 1)

# add header to matrix results file
echo -e 'Sample\t'$contingList > $resultsFile

# loop over each sample
for f in $inputsPath"/"*"_trimmed.vcf"; do
	# retrieve sample tag
	sampleTag=$(basename $f1 | sed "s/\_trimmed\.vcf$//g")
	# add sample tag to matrix row
	echo $sampleTag >> $resultsFile
	# loop over each marker
	for i in $contingList; do
		# create file with current marker variants
		cat $f | grep $i > $f"."$i".tmp.txt"
		# initialize alleles
		firstAlleles="NULL"
		secondAlleles="NULL"
		# loop over each line in the sample vcf
		while read line; do
			# retrieve the CHROM
			chrom=$(echo $line | cut -f 1)
			# retrieve the POS
			pos=$(echo $line | cut -f 2)
			# retrieve the REF
			ref=$(echo $line | cut -f 4)
			# retrieve the ALT
			alt=$(echo $line | cut -f 5)
			# retrieve the GT alleles and translate encodings
			firstGT=$(echo $line | cut -f 10 | cut -d ":" -f 1 | cut -d "/" -f 1 | sed "s/0/$ref/g" | sed "s/1/$alt/g")
			secondGT=$(echo $line | cut -f 10 | cut -d ":" -f 1 | cut -d "/" -f 2 | sed "s/0/$ref/g" | sed "s/1/$alt/g")
			# add POS to GT alleles
			firstGT=$firstGT"("$pos")"
			secondGT=$secondGT"("$pos")"
			# append GT to lists of allele variants
			firstAlleles=$firstAlleles","$firstGT
			secondAlleles=$secondAlleles","$secondGT
		done < $f"."$i".tmp.txt"
		# output allele lists to the results matrix
		echo -en $firstAlleles'\t'$secondAlleles >> $resultsFile
	done
done

# status message
echo "Analysis conplete!"

#!/bin/bash

# script to format input data for the SSR pipelines
# usage: bash inputs_formatting.sh
# usage Ex: qsub inputs_formatting.sh

# load software module
module load bio

# retrieve reference absolute path for alignment
ref=$(grep "reference:" "inputs_ssr_pipeline.txt" | tr -d " " | sed "s/reference://g")

# status message
echo "Formatting started..."

# first, make sure the ssr info and primer sequence files have been converted from csv to txt
bash ssr_info_formatting.sh
bash ssr_regions_formatting.sh
bash primer_info_formatting.sh

# make sure the reference has been indexed
bwa index $ref

# status message
echo "Formatting complete!"

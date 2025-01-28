# GBCF Bioinformatics Analysis Report

- Project: JRS Chestnut SSR Analysis
- Contact: Jeanne Romero-Severson
- Analyst: Elizabeth Brooks
- Date: 28 January 2025

## Methods

Custom BASH and python scripts used in the analysis are available on GitHub (https://github.com/ElizabethBrooks/SSRToolDrawer). Raw sequences were trimmed of adapters with Trimmomatic version 0.39 (Bolger et al., 2014) and assessed for quality with FastQC v0.11.8 (Andrews, 2010). Trimmed sequences were aligned to the reference marker contigs using the BWA software package version 0.7.17-r1188 (Li et al., 2009). Corresponding alignments were sorted and sample read groups applied with SAMtools and BCFtools versions 1.9 (Danecek et al., 2021). Primer sequences were soft masked using BAMClipper v1.0.0 (Au et a., 2017). Variants were called and filtered using BCFtools, then trimmed of SSR regions using BEDTools v2.30.0 (Quinlan & Hall 2010). Custom Python2 scripts created previously by Joseph Sarro were used to measure and filter SSR lengths. Next, a formatted matrix of SSR lengths was created for the basic workflow. A custom BASH script was created by Elizabeth Brooks to format the matrix of variants resulting from the SNP calling workflow. The entire workflow of analysis was combined into a pipeline using custom BASH scripts created by Elizabeth Brooks.

## References

1. Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data[Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
2. Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for IlluminaSequence Data. Bioinformatics, btu170
3. Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168] Twelve years of SAMtools and BCFtools Petr Danecek, James K Bonfield, Jennifer Liddle, John
4. Marshall, Valeriu Ohan, Martin OPollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li GigaScience, Volume 10, Issue 2, February 2021, giab008,https://doi.org/10.1093/gigascience/giab008
5. Au, C., Ho, D., Kwong, A. et al. BAMClipper: removing primers from alignments to minimize false-negative mutations in amplicon next-generation sequencing. Sci Rep 7, 1567 (2017). https://doi.org/10.1038/s41598-017-01703-6
6. Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics (Oxford, England), 26(6), 841–842. https://doi.org/10.1093/bioinformatics/btq033

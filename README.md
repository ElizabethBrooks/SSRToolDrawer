# SSRTools

## Installation

<b>First,</b> download this GitHub repository using the git clone command in the terminal as follows.

To download the code onto a local computer or server space, click the green <i>< > Code</i> button and copy the link. Then, using the HTTPS web URL in the terminal:

git clone https://github.com/ElizabethBrooks/SSRToolDrawer.git

Alternatively, using SSH:

git clone git@github.com:ElizabethBrooks/SSRToolDrawer.git

The latest version of this application may also be downloaded from this repository by clicking the green <i>< > Code</i> button near the top of the page, and then clicking <i>Download ZIP</i>.

## Running

1. Update the **outputs** path in the *inputs_sst_pipeline.txt* file. Usually the final directory in the outputs path is updated to the current month and year like so: SSRAnalysis\_**MonthYear** (e.g., SSRAnalysis_May2024).
2. Add a file to the InputData directory that contains the **run**, **ID**, and **pairedReads**. The file should be named with the run number like so: inputs_run**Number**.txt (e.g., inputs_run1.txt). For example input data file formatting, refer to the *inputs_run1.txt* file.
3. To run the "Basic" workflow, move to the **Basic** directory and run the *ssr_pipeline_basic.sh* script by following the usage notes at the top of the file (e.g., qsub ssr_pipeline_basic.sh inputs_run8.txt).
4. To run the "SNP Calling" workflow, move to the **SNP_Calling** directory and run the *ssr_pipeline_SNP.sh* script by following the usage notes at the top of the file (e.g., qsub ssr_pipeline_SNP.sh run1 run2 run3 run4 run5 run6 run7 run8).

### Note

The *ssr_pipeline_basic.sh* (lines 77 and 78) file in the Basic directory and *ssr_pipeline_subset_SNP.sh* (lines 64 and 65) file in the Scripts subdirectory of the SNP_Calling directory uses the bashrc and python2 environment setup by the ebrooks5 user.

## Workflows

There are two possible workflows.

### Basic
The sequence capture data should first be assessed for quality with fastqc, trimmed with trimmomatic, and mapped with bwa. After mapping, a data matrix consisting of SSR information can be generated using the following set of python2 scripts and workflow.

- \<QC with fastqc, trim with trimmomatic, and map with bwa\>
- python2 GapGenes.v3.py -sam \<SAM file\> -C \<SSR info file\> -P/S
- python2 SnipMatrix.py \<Output file from GapGenes\>
- python2 Format_Matrix.py

### SNP calling
The sequence capture data should first be assessed for quality with fastqc, trimmed with trimmomatic, and mapped with bwa. A python2 script can then be run to print out sequences that flank SSRs at least 50bp in both directions. Next, the aligned data can be filtered to only include reads mapped within 50bp of a SSR before SNPs are called and VCF files generated with bcftools. Finally, a data matrix consisting of SSR information can be generated using the following set of python2 scripts and workflow.
    
- \<QC with fastqc, trim with trimmomatic, and map with bwa\>
- python2 SamIAm.py -sam \<SAM file\> -C \<SSR info file\> -p \<yes/no\>
- \<Filter to only include reads mapped within 50bp of a SSR\>
- \<Call SNPs with bcftools\>
- python2 Format_VCF-Matrix.py


## Description of Tools
The following scripts were generated by [Joe Sarro](https://github.com/jsarro13/SSRTools) to form a pipeline that takes sequence capture data and generates a data matrix consisting of Short Sequence Repeat (SSR) information. 

In addition to the novel scripts, two [previously written scripts](https://git.io/JTYMb) were used. FindReadLengths.py was used to verify our sequence capture had worked correctly, and SnipVCF.py to filter downstream SNPs. 

### sardinesandwich.py 
This script can be used as a quality control measure after sequence alignment. It will output four files. One file containing reads that flank 50bp or more both upstream and downstream of the SSR, a file containing reads that only flank 50bp or more upstream of the SSR, a file containing reads that flank 50bp or more only downstream of the SSR, and reads that do neither. Files contain three columns, the read length, flanking length to the left of the SSR, and flanking length to the right of the SSR. The script will also print all counts to terminal.

### SamIAm.py
This script will print out sequences that flank and SSR at least 50bp in both directions. These sequences can be used for downstream SNP calling 

### GapGenes.py 
Measures the length of an SSR after sequence alignment. Two output files are generated. The first prints the SSR, repeat number, and both flanking ends. A second file contains a matrix of count information of SSR lengths.

### SNPMatrix.py
This script will filter the matrix output of GapGenes.py. It filters with the following criteria: 1) A contig can contain no more than 3 SSRs 2) An SSR must have at least 4 reads. 3) SSR two and three may not be less than 25% the most abundant SSR.

### FormatMatrix.py
Formats the output of SNPMatrix.py into a format that can be read by JoinMap.


## Workflow Input Files

### SSR Information
- /afs/crc.nd.edu/group/genomics/lando/JRS-chestnut_02-2022_miseq/SSR_info.txt

### SSR Sequences
- /afs/crc.nd.edu/group/genomics/lando/JRS-chestnut_02-2022_miseq/finalset-2.fa

### Sequence Capture Data
- run1: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/190920_Chestnut_Multiplex-PCR_Plate02
- run2: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/211111_JB-688_Romero-Severson-Chestnut
- run3: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/220128_JB-861_Romero-Severson-Chestnut/Fastq
- run4: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/220408_JB-998_ChestnutSSR_Round4/Fastq
- run5: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/220706_ChestnutSSR_Round5
- run6: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/230207-230130_JRS-ChestnutSSR-Set6
- run7: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/230719_ChestnutSSR_Round7
- run8: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/240112_DS-2529_ChestnutSSR_Round8

#### Note

The plate 6 sets were split across two runs:
- /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/230130_JRS-ChestnutSSR-Set6
- /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/230207_JRS-ChestnutSSR-Set6

## Past Analysis Data Sets

Original test: 
- /afs/crc.nd.edu/group/genomics/yoda/JRS_MiSeq_7-27-17

1st chestnut:
- /afs/crc.nd.edu/group/genomics/Leia/JRS_Chestnut_11_2019
- /afs/crc.nd.edu/group/genomics/Leia/JRS_Chestnut_11_2019-2

Walnut:
- /afs/crc.nd.edu/group/genomics/Leia/JRS_Black_Walnut_03_2020-1
- /afs/crc.nd.edu/group/genomics/Leia/JRS_Black_Walnut_03_2020-2

Last chestnut:
- /afs/crc.nd.edu/group/genomics/lando/JRS-chestnut_2021_miseq
- /afs/crc.nd.edu/group/genomics/lando/JRS-chestnut_02-2022_miseq

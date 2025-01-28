### Analysis Pipeline

#### Prep

For each run:

1. QC (FastQC)
2. Trimming (Trimmomatic)
3. Mapping (BWA)

#### SSR Lengths

For each run:

1. Measuring SSR lengths (GapGenes.v3.py)
2. Filtering the matrix of SSR lengths (SnipMatrix.py)
3. Formatting of SSR lengths matrix for JoinMap (Format_Matrix.py)

#### SNP Calling

For each run:

1. Retrieval of sequences flanking SSRs at least 50bp in both directions (SamIAm.py)
2. Sorting and removal of pcr duplicates (SAMtools)
3. Filtering of sequences to keep only unique read alignments (SAMtools)
4. Clipping to soft mask primer sequences (BAMClipper)
5. Addition of sample read groups (SAMtools)

Across all runs:

6. Variant calling (BCFtools)
7. Variant filtering (BCFtools)
8. Variant trimming to remove SSR regions (BEDTools)
9. Formatting of variants matrix (variantMatrix_bcftools.sh)

#### Notes

1. The basic SSR lengths workflow is performed on a per-run basis
2. The SNP calling workflow is performed with the samples from all of the available runs
(run1 to run9)
3. The run number for each sample is appended to the end of each sample name (e.g.,
SAMPLE1_run1) in the results matrix from the SNP calling workflow


### Workflow Input Files

#### SSR Information
- /afs/crc.nd.edu/group/genomics/lando/JRS-chestnut_02-2022_miseq/SSR_info.txt

#### SSR Sequences
- /afs/crc.nd.edu/group/genomics/lando/JRS-chestnut_02-2022_miseq/finalset-2.fa

#### Sequence Capture Data - Working
- run1: /afs/crc.nd.edu/group/genomics/Mando/ChestnutSSR_Round1
- run2: /afs/crc.nd.edu/group/genomics/Mando/ChestnutSSR_Round2
- run3: /afs/crc.nd.edu/group/genomics/Mando/ChestnutSSR_Round3
- run4: /afs/crc.nd.edu/group/genomics/Mando/ChestnutSSR_Round4
- run5: /afs/crc.nd.edu/group/genomics/Mando/ChestnutSSR_Round5
- run6: /afs/crc.nd.edu/group/genomics/Mando/ChestnutSSR_Round6
- run7: /afs/crc.nd.edu/group/genomics/Mando/ChestnutSSR_Round7
- run8: /afs/crc.nd.edu/group/genomics/Mando/ChestnutSSR_Round8
- run9: /afs/crc.nd.edu/group/genomics/Mando/ChestnutSSR_Round9
- run10: /afs/crc.nd.edu/group/genomics/Mando/ChestnutSSR_Round10

#### Sequence Capture Data - Original
- run1: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/190920_Chestnut_Multiplex-PCR_Plate02
- run2: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/211111_JB-688_Romero-Severson-Chestnut
- run3: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/220128_JB-861_Romero-Severson-Chestnut
- run4: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/220408_JB-998_ChestnutSSR_Round4
- run5: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/220706_ChestnutSSR_Round5
- run6: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/230207-230130_JRS-ChestnutSSR-Set6
- run7: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/230719_ChestnutSSR_Round7
- run8: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/240112_DS-2529_ChestnutSSR_Round8
- run9: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/Romero-Severson_DS-2862_ChestnutSSR_Round9_240612_M00314_LL7GG
- run10: /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/Romero-Severson_DS-3011_240824_M00314

##### Note

The plate 6 sets were split across two runs:
- /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/230130_JRS-ChestnutSSR-Set6
- /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/230207_JRS-ChestnutSSR-Set6

Combined plate 6 sets:
- /afs/crc.nd.edu/group/genomics/DEATHSTAR/MiSeq/230207-230130_JRS-ChestnutSSR-Set6

### Past Analysis Data Sets

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

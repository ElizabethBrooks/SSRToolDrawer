#!/usr/bin/python
#Fromat Matrix  03-2020
#Written by Joseph Sarro
#Format read 1 matrix file to  
from __future__ import division
import numpy as np
import os
import argparse
from operator import itemgetter
# contig (marker) ID list
contig_list=["001-QR_002", "002-CD_001", "003-CD_003", "004-CD_005", "005-CD_011", "006-CD_013", "007-CD_014", "008-CD_015", "009-CD_017", "010-CD_019", "011-CD_020", "012-CD_024", "013-CD_026", "014-CD_028", "015-CD_033", "016-CD_042", "017-CD_043", "018-CD_049", "019-CD_053", "020-CD_056", "021-CD_057", "022-CD_058", "023-CD_059", "024-CD_061", "025-CD_063", "026-CD_066", "027-CD_067", "029-CD_077", "030-CD_082", "031-CD_084", "032-CD_087", "033-CD_088", "035-CD_095", "036-CD_100", "037-CD_101", "038-CD_104", "039-CD_105", "047-CD_122", "049-QR_003", "050-CM_001", "051-CM_003", "052-CM_007", "053-CM_014", "057-QR_012", "058-QR_041", "059-QR_046", "060-QR_055", "061-QR_056", "062-QR_059", "063-QR_062", "064-QR_069", "065-QR_077", "066-QR_088", "067-QR_092", "068-QR_098", "069-QR_122", "070-QR_150", "071-QR_181", "072-QR_186", "073-QR_192", "074-QA_001", "075-QA_004"]
# sample ID list
sample_list=["FIND_ME_REPLACE_ME"]
# set output matrix file
G= open("VCF_Matrix.txt",'w')
print >> G, "Sample	",
for i in contig_list:
	# add contig ID to the matrix header
	print >> G, i,"		",  
print >> G
# loop over each sample in the list
for i in sample_list:
	# add sample ID to the first column of the matrix
	print >> G, i, 
	# loop over each contig in the list
	for j in contig_list:
		# flag to check if a cell has been made for the current contig
		found= False
		# flag to check the previous contig in the loop
		prev="NA"
		# set sample vcf file
		ta = open(i+".noHeader.vcf", 'r')
		for lines in ta:
			fields = lines.split("	")
			contig=fields[0]
			# the reference position, with the 1st base having position 1
			pos=fields[1]
			# comma separated list of alternate non-reference alleles
                        alt=fields[4]
            # phred-scaled quality score
                        qual=fields[5]
            # check if the current contig has been processed and is in the contig list
			if prev.rstrip(" ")==contig.rstrip(" ") and j.rstrip(" ")==contig.rstrip(" "):
				#print "here"
				if float(qual)>=30:
					# add alt and pos to the current cell in the sample row
                                        print >> G,",",alt,"(",pos,")",
            # check if the current contig is in the contig list
			elif j.rstrip(" ")==contig.rstrip(" "):
				prev=contig
				if float(qual)>=30:
					found=True
					# add alt and pos to a new cell in the sample row
                                        print >> G,"	",alt,"(",pos,")",
		# check if a new cell has been made for the current contig
		if found == False:
			# output NULL to the new cell in the sample row
			print >> G,"	NULL",  	
		ta.close()
	print >> G
# close and tab format output matrix
G.close()
os.system("sed -i -e  's/\t/_NULL_/g' VCF_Matrix.txt")
os.system("sed -i -e  's/\s//g' VCF_Matrix.txt")
os.system("sed -i -e  's/_NULL_/\t/g' VCF_Matrix.txt")

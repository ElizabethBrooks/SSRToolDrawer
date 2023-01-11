#!/usr/bin/python
#Fromat Matrix  03-2020
#Written by Joseph Sarro
#Format read 1 matrix file to  
from __future__ import division
import numpy as np
import os
import argparse
from operator import itemgetter
contig_list=["001-QR_002", "002-CD_001", "003-CD_003", "004-CD_005", "005-CD_011", "006-CD_013", "007-CD_014", "008-CD_015", "009-CD_017", "010-CD_019", "011-CD_020", "012-CD_024", "013-CD_026", "014-CD_028", "015-CD_033", "016-CD_042", "017-CD_043", "018-CD_049", "019-CD_053", "020-CD_056", "021-CD_057", "022-CD_058", "023-CD_059", "024-CD_061", "025-CD_063", "026-CD_066", "027-CD_067", "029-CD_077", "030-CD_082", "031-CD_084", "032-CD_087", "033-CD_088", "035-CD_095", "036-CD_100", "037-CD_101", "038-CD_104", "039-CD_105", "047-CD_122", "049-QR_003", "050-CM_001", "051-CM_003", "052-CM_007", "053-CM_014", "057-QR_012", "058-QR_041", "059-QR_046", "060-QR_055", "061-QR_056", "062-QR_059", "063-QR_062", "064-QR_069", "065-QR_077", "066-QR_088", "067-QR_092", "068-QR_098", "069-QR_122", "070-QR_150", "071-QR_181", "072-QR_186", "073-QR_192", "074-QA_001", "075-QA_004"]
sample_list=["FIND_ME_REPLACE_ME"]
G= open("VCF_Matrix.txt",'w')
print >> G, "Sample	",
for i in contig_list:
	print >> G, i,"		",  
print >> G
for i in sample_list:
	#os.system("cp "+i+".sam.Matrix.1.txt.trimmed.txt temp")
	#os.system("sed temp -e '1,1d' > tempA")
	#os.system("sed -i -e 's/([0-9]\+)//g' tempA")
	#os.system("rm temp")
	##ta = open("tempA", 'r')
	print >> G, i, 
	for j in contig_list:
		found= False
		prev="NA"
		ta = open(i+".noHeader.vcf", 'r')
		#print ta
		for lines in ta:
			fields = lines.split("	")
			contig=fields[0]
		#	print contig
			pos=fields[1]
                        alt=fields[4]
                        qual=fields[5]
			if prev.rstrip(" ")==contig.rstrip(" ") and j.rstrip(" ")==contig.rstrip(" "):
				#print "here"
				if float(qual)>=30:
                                        print >> G,",",alt,"(",pos,")",
			elif j.rstrip(" ")==contig.rstrip(" "):
				prev=contig
                        	#ssr = ssr.rstrip("\n")
                       		#ssrs=ssr.split(',')
				#print >> G
				if float(qual)>=30:
					found=True
                                        print >> G,"	",alt,"(",pos,")",
				#else:
				#	print  >> G, "	",ssrs[0],"	",ssrs[1],"	",ssrs[2],
		if found == False:
			print >> G,"	NULL",  	
		ta.close()
	print >> G
#print >> G, "Total_Library","  ",total_count_A,"       ",matched_count_A,"     ",'{:.1%}'.format(matched_count_A/total_count_A)
#os.system("rm tempA")
G.close()
os.system("sed -i -e  's/\t/_NULL_/g' VCF_Matrix.txt")
os.system("sed -i -e  's/\s//g' VCF_Matrix.txt")
os.system("sed -i -e  's/_NULL_/\t/g' VCF_Matrix.txt")

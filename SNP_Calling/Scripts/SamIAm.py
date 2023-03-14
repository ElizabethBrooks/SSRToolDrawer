#!/usr/bin/python
#SamIAm 02-2018
#Written by Joseph Sarro
# Edited by Elizabeth Brooks
# Edited Feb 2023
#This script inputs a sam file and a list of coordinates of SSRs and prints reads that span the SRR 50 bp in each direction.
import sys
import argparse
#Check to see if reads are single end or paired end.  Tye 2 as an argument for paired end and 1 for single end, default is 1.
SP = 1
#files to be opened and written.  F needs to be changed to match your SAM file before running as does datafile below.  
samplename =' '
parser = argparse.ArgumentParser(description='Filters a SAM file with reads that span an SSR at least 50bp in each direction.')
parser.add_argument('-sam', type=str,help='Enter a SAM file')
parser.add_argument('-C', type=str, required=True, help='Enter a coordinates file')
parser.add_argument('-p', dest='pairs',type=str,default='no',help='Enter yes if these are paired reads (Defualt=no)')
args = parser.parse_args()
if args.pairs == 'yes':
	SP=2
print SP
samplename= args.sam
f = open(samplename, 'r')
#f = open('S6b.sam', 'r')
g = open(samplename+'.filter50.sam','w')
seqHitInfo = open(samplename+'.hitInfo','w')
count=0
mapcount=0
contigcount=0
hitcount=0
hitrightcount=0
hitleftcount=0
hitnonecount=0
whatcount=0
#Are the reads aired for single end?
if SP ==2:
	linecount=1
	readlen=0
	start=0
	for lines in f:
		fields = lines.split("	")
		sq=fields[0]
		#ignore lines with headers
		if sq != "@SQ" and sq != "@RG" and sq != "@PG":
			linecount+=1
			if linecount % 2 != 0:
				count +=1
				start=fields[3]
				field1= fields[9]
				readlen=len(field1)
			else:
                                endstart=fields[3]
                                field1= fields[9]
                                endreadlen=len(field1)
				end = int(endstart)+endreadlen-1
				contig=fields[2]
				#all mapped reads
				if contig !='*':
					mapcount+=1
					#A file containing adresses of all SSRs, datafile needs to be changed to your file name.  Will mark as a hit if the read spands 50 bp to the left and right of the SRR.  
					#Files containing the length of the sequence, number of bp spanning left of the SSR, and right of the SSR, will be printed to appropriate output files. 
					datafile = file(args.C)
					for line in datafile:
						if contig in line:
							contigcount+=1
							hit=line.split("	")
							hitstart=hit[1]
							hitend=hit[2]
							left5span= int(hitstart)-int(start)
							right3span=int(end)-int(hitend)
							if int(start) <= (int(hitstart)-50) and int(end) >= (int(hitend)+50):
								hitcount +=1
								print >> g, lines,
								print >> seqHitInfo, readlen, left5span, right3span
		elif sq == "@SQ" and sq == "@RG" and sq == "@PG":
			print >> g, lines,	
else:
#for all lines in the sam file
	for lines in f:
		fields = lines.split("	")
		sq=fields[0]
		#ignore lines with headers
		if sq != "@SQ" and sq != "@RG" and sq != "@PG":
			count +=1
			start=fields[3]
			field1= fields[9]
			readlen=len(field1)
			end = int(start)+readlen-1
			contig=fields[2]
			#all mapped reads
			if contig !='*':
				mapcount+=1
				#A file containing adresses of all SSRs, datafile needs to be changed to your file name.  Will mark as a hit if the read spands 50 bp to the left and right of the SRR.  
				#Files containing the length of the sequence, number of bp spanning left of the SSR, and right of the SSR, will be printed to appropriate output files. 
				datafile = file(args.C)
				for line in datafile:
					if contig in line:
						contigcount+=1
						hit=line.split("	")
						hitstart=hit[1]
						hitend=hit[2]
						# calculate the number of bases spanning left (start) and right (end)
						left5span= int(hitstart)-int(start)
						right3span=int(end)-int(hitend)
						# hit if there are at least 50 bases spanning to the left (start) and right (end) of a SSR
						if int(start) <= (int(hitstart)-50) and int(end) >= (int(hitend)+50):
							# increment hit count
							hitcount +=1
							# output sequence hit
							print >> g, lines,
							# output sequence hit info
							print >> seqHitInfo, readlen, left5span, right3span
		else:
			print >> g, lines,

# print total hit number
print "Total sequences with 50 bp spanning to the left and right of a SRR: ", hitcount

#######################################################################
##
## FilterSequences.py
##
## Version 2.0 -- 6 April 2013
##
## Created by Michael Sorenson
## Copyright (c) 2011-2013 Boston University. All rights reserved.
##
## Dependencies: usearch (version 5 in $PATH as "usearch")
##
## This program is written for execution in the python (version 3)
## language. It is free and distributed WITHOUT warranty; without
## even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.
##
#######################################################################

import sys
import os, subprocess, heapq

#qseq_filename = "2270-37420_ER_TNC_S1_L002_R1_001tempsort.qseq"
#qseq_file1 = open(qseq_filename,'r')

for qseq_filename in glob.glob('*tempsort.qseq'):
	print(qseq_filename)

	#sort by quality
	print('Sorting qseq file by quality')
	command=('sort -k13,13nr -k1,1d -k4,4n -k5,5n -k6,6n '+qseq_filename+' -o '+qseq_filename.replace('.qseq','sort.qseq'))
	p = subprocess.Popen(command, shell=True)
	sts = p.wait()

	#make fasta file
	print('Making fasta file for uClust')
	command = ("awk '//{print(\">\"$1\":\"$4\":\"$5\":\"$6\"\\n\"$9)}' "+qseq_filename.replace('.qseq','sort.qseq')+" > "+qseq_filename.replace('.qseq','.fasta'))
	print(command)
	p = subprocess.Popen(command, shell=True)
	sts = p.wait()

	#run uClust on fasta file
	print('Running uclust on fasta file')
	command = ('/Volumes/VIREO/SNP_Test/usearch -cluster_fast ' + qseq_filename.replace('.qseq','.fasta') + ' -id 0.9 -centroids nr.fasta -uc ' + qseq_filename.replace('.qseq', '.uc'))
	p = subprocess.Popen(command, shell=True)
	sts = p.wait()

	##get cluster depths
	print('Getting cluster depths')
	command=("awk '/C\t/{print$2,$3}' "+qseq_filename.replace('.qseq','.uc')+" > "+qseq_filename.replace('.qseq','.clusterdepth'))
	p = subprocess.Popen(command, shell=True)
	sts = p.wait()

	depths=[]
	clusterdepth_file = open(qseq_filename.replace('.qseq','.clusterdepth'),'r')
	for line in clusterdepth_file:
	clusterdata = line.split()
	#print("cluster data is: ")
	#print(clusterdata)
	depths.append(int(clusterdata[1]))
	#print(depths)


	#get cluster numbers and add to sorted file
	ucresults_file = open(qseq_filename.replace('.qseq','.uc'),'r')
	for i in range(8):
	ucresults_file.readline()

	qseq_file1 = open(qseq_filename.replace('.qseq','sort.qseq'), 'r')
	qseq_file2 = open(qseq_filename.replace('.qseq','temp.qseq'), 'w')

	for line in ucresults_file:
	data = line.split()
	if data[0] == "C":
		break
	else:
		line2 = qseq_file1.readline()
		qseq_file2.write(line2[:line2.rfind('\n')]+'\t'+data[1]+'\n')
	qseq_file2.close()

	#sort by cluster number and quality
	print('Sorting qseq file by quality')
	command=('sort -k15,15n -k13,13nr -k1,1d -k4,4n -k5,5n -k6,6n '+qseq_filename.replace('.qseq','temp.qseq')+' -o '+qseq_filename.replace('.qseq','sort2.qseq'))
	p = subprocess.Popen(command, shell=True)
	sts = p.wait()

	#evaluate each cluster and filter seqs
	qseq_file1 = open(qseq_filename.replace('.qseq','sort2.qseq'), 'r')
	outf1 = open(qseq_filename.replace('.qseq','X.qseq'),'w')
	outf2 = open(qseq_filename.replace('.qseq','F.qseq'),'w')

	line_counts = [sum(depths[:i+1]) for i in range(len(depths))]
	lines = qseq_file1.readlines()
	#print(type(lines))


	filtertest = [0]*10
	previous = 0
	print("Looping through")
	for i in range(len(line_counts)):
	cluster_data = []
	for j in range(previous,line_counts[i]):
		cluster_data.append(lines[j])
		#print(lines[j])
	#for j in cluster_data:
		#print(j)
	#print(cluster_data)

	temp = []        
	for k in cluster_data:
		k = k.split()
		temp.append([int(k[11]),float(k[12])])
	count_quals = list(zip(*temp))
	print("Count quals:")
	print(type(count_quals[0]))
	print(count_quals[0])

	while(count_quals[0] ==1):
		if counts_quals[1,0] < 20 or cluster_data[0,8].find('NNNNNNNNNN') > 0:
			filtertest[7]+=1
			outf1.write('\t'.join(cluster_data[0])+'\n')
		elif (cluster_data[0][13]) > 10 and(cluster_data[0][12]) < 30:
			filtertest[8]+=1
			outf1.write('\t'.join(cluster_data[0])+'\n')
		else:
			filtertest[9]+=1
			outf2.write('\t'.join(cluster_data[0][:13])+'\n')

	while(count_quals[0] != 1):
		trim = False
		top_3 = heapq.nlargest(3,count_quals[0])
		#print("Top 3:")
		#print(top_3[0])
		#print(type(count_quals[0]))
		#print(type(top_3[0]))
	
		if (top_3[0]/(sum(count_quals[0][0:]))) > 0.80:
			if ((top_3[0] > 29 and top_3[1])< 3):
				filtertest[0]+=1
				for k in cluster_data:
					if int(k[11]) < 3:
						outf1.write('\t'.join(k)+'\n')
					else:
						outf2.write('\t'.join(k[:13])+'\n')
			else:
				filtertest[1]+=1
				for k in cluster_data:
					if float(k[12])<30:
						outf1.write('\t'.join(k)+'\n')
					else:
						outf2.write('\t'.join(k[:13])+'\n')
		if sum(top_3[:2])/sum(count_quals[0][0:]) > 0.80 and len(top_3) > 2:
			if sum(top_3[:2]) > 29 and top_3[2] < 3 and top_3[0]/top_3[1] < 1.6:
				filtertest[2]+=1
				for k in cluster_data:
					if int(k[11]) < 3:
						outf1.write('\t'.join(k)+'\n')
					else:
						outf2.write('\t'.join(k[:13])+'\n')
			else:
				filtertest[3]+=1
				for k in cluster_data:
					if float(k[12])<30:
						outf1.write('\t'.join(k)+'\n')
					else:
						outf2.write('\t'.join(k[:13])+'\n')
		while sum(count_quals[0]) > 10 and max(count_quals[1]) > 35:
			filtertest[4]+=1
			for k in cluster_data:
				if float(k[12])<20:
					outf1.write('\t'.join(k)+'\n')
				else:
					outf2.write('\t'.join(k[:13])+'\n')
		while sum(count_quals[0]) > 10 and max(count_quals[1]) > 30:
			filtertest[5]+=1
			for k in cluster_data:
				if float(k[12])<10:
					outf1.write('\t'.join(k)+'\n')
				else:
					outf2.write('\t'.join(k[:13])+'\n')
		else:
			filtertest[6]+=1
			for k in cluster_data:
				outf2.write('\t'.join(k[:13])+'\n')

					 
	previous = line_counts[i]
	
	print(filtertest)

	outf1.close()
	outf2.close()

	##cleanup
	command=('rm '+qseq_filename.replace('.qseq','temp.qseq')+' '+qseq_filename.replace('.qseq','sort2.qseq')+' '+qseq_filename.replace('.qseq','sort.qseq')+' '+qseq_filename.replace('.qseq','.clusterdepth')+' '+qseq_filename.replace('.qseq','.uc90')+' '+qseq_filename.replace('.qseq','.fas'))
	p = subprocess.Popen(command, shell=True)
	sts = p.wait()
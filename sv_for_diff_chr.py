'''
File: sv_for_diff_chr.py
Author: Rendong
Description: detect inter-chromosomal translocation supported by pairs mapped to different chromosome
History: Thu Feb 23 15:57:35 EST 2012
'''
import os,sys
import re

def intersect(my_range,start):
	"""if overlap between two ranges, update my_range and return True"""
	if start < my_range[1]:
		if start+read_len+search_len < my_range[1]:
			return False
		else:
			my_range[1] = start
			return True
	elif start > my_range[2]+search_len:
		return False
	else:
		my_range[2] = start+read_len
		return True

if len(sys.argv)<2:
	print "Usage:python sv_for_diff_chr.py diffchr_pairs.txt output.txt read_len(optional, default 100) search_len(optional, default 1000) cluster_cutoff(optional, default 1)"
	sys.exit(0)
elif len(sys.argv) == 3:
	read_len = 100
	search_len = 1000
	cut_off = 1
else:
	read_len = int(sys.argv[3])
	search_len = int(sys.argv[4])
	cut_off = len(sys.argv[5])

left_range = ['chr0',0,1e20]
right_range = ['chr0',0,1e20]
pair_set = set() 
cluster = []
ifp = open(sys.argv[1])
ofp = open(sys.argv[1]+'.sv.tmp.txt','w')
for i,line in enumerate(ifp):
	items = line.rstrip().split()
	if items[0] in pair_set:
		continue
	else:
		pair_set.add(items[0])
	#get cluster for pairs
	if items[2] == left_range[0] and (int(items[3]) >= left_range[1] and int(items[3]) <= left_range[2]+search_len):	#left end match
		if items[6] == right_range[0] and intersect(right_range,int(items[7])):	#right end match
			cluster.append(i+1)
			left_range[2] = int(items[3])+read_len
			continue
	if len(cluster) >= cut_off:
		print >> ofp,'inter-chromosome\t'+str(len(cluster))+'\t'+','.join(map(str,cluster))+'\t'+'\t'.join(map(str,left_range))+'\t'+'\t'.join(map(str,right_range))
	cluster = [i+1]
	left_range = [items[2],int(items[3]),int(items[3])+read_len]
	right_range = [items[6],int(items[7]),int(items[7])+read_len]

if len(cluster) >= cut_off:
	print >> ofp,'inter-chromosome\t'+str(len(cluster))+'\t'+','.join(map(str,cluster))+'\t'+'\t'.join(map(str,left_range))+'\t'+'\t'.join(map(str,right_range))
ifp.close()
ofp.close()

# adding genomic annotation information to detected sv file.
ifp = open(sys.argv[1]+'.sv.tmp.txt')
ofp = open(sys.argv[2],'w')
print >> ofp, 'type\tpair_counts\trelated_pairs_ID\tleft_chr\tleft_start\tleft_end\tright_chr\tright_start\tright_end\tleft_ann\tright_ann\timportant'
ofp1 = open(sys.argv[1]+'.left.bed','w')
ofp2 = open(sys.argv[1]+'.right.bed','w')
for line in ifp:
	items = line.rstrip().split()
	p = re.match('chr(.*)',items[3])
	if p.group(1) == 'X':
		items[3] = '23'
	elif p.group(1) == 'Y':
		items[3] = '24'
	else:
		items[3] = p.group(1)
	p = re.match('chr(.*)',items[6])
	if p.group(1) == 'X':
		items[6] = '23'
	elif p.group(1) == 'Y':
		items[6] = '24'
	else:
		items[6] = p.group(1)
	print >> ofp1,'\t'.join(items[3:6])
	print >> ofp2,'\t'.join(items[6:9])
ifp.close()
ofp1.close()
ofp2.close()
os.system('perl addcatann.pl human 10 2 2 '+sys.argv[1]+'.left.bed 0 hg19.refseq.unique.annot '+sys.argv[1]+'.left.out')
os.system('perl addcatann.pl human 10 2 2 '+sys.argv[1]+'.right.bed 0 hg19.refseq.unique.annot '+sys.argv[1]+'.right.out')
ifp = open(sys.argv[1]+'.sv.tmp.txt')
ifp1 = open(sys.argv[1]+'.left.out')
ifp2 = open(sys.argv[1]+'.right.out')
while True:
	line = ifp.readline().rstrip()
	if not line:
		break
	line1 = ifp1.readline().rstrip().split('\t')
	line2 = ifp2.readline().rstrip().split('\t')
	l_ann = ','.join(line1[3:])
	r_ann = ','.join(line2[3:])
	if 'intergenic' in l_ann or 'intergenic' in r_ann:
		important = 0
	else:
		important = 1
	print >> ofp, line+'\t'+l_ann+'\t'+r_ann+'\t'+str(important)
ifp.close()
ifp1.close()
ifp2.close()
ofp.close()
os.system('rm '+sys.argv[1]+'.sv.tmp.txt')
os.system('rm '+sys.argv[1]+'.left*')
os.system('rm '+sys.argv[1]+'.right*')
print 'finish clustering pairs'

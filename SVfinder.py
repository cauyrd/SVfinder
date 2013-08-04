'''
File: SVfinder.py
Author: Rendong Yang
Email: cauyrd@gmail.com
Description: detect Structural Variations supported by discordant read pairs
History: 1.Thu Feb 23 15:57:35 EST 2012
		 2.Wed Jul 31 16:56:57 EDT 2013
'''
import os,sys
import re,getopt
import HTSeq

def extend_interval(iv1,iv2):
	"""newest HTSeq package has extend_to_include function do the same work"""
	if iv1.strand != iv2.strand or iv1.chrom != iv2.chrom:
		print "Genomic intervals with different chrom or strand cannot be merged!"
		sys.exit(1)
	iv1.start = min(iv1.start, iv2.start)
	iv1.end = max(iv1.end, iv2.end)

def usage():
	"""showing help information"""
	print 'Usage:'
	print 'python SVfinder.py -i <discordant_pairs.sam> -o <output_file> [opts]'
	print 'Opts:'
	print ' -n <int>	:cutoff of number of discordant pais to define a cluster  (default:2)'
	print ' -l <int>	:extention length to join overlaped reads together (default:1000)'
	print ' -r <int>	:read length (default:100)'
	print ' -g <int>	:gene annotation file (default:hg19.ucsc.gene.txt)'
	print ' -h      	:produce this menu'

# parameters parsing.
lread = 100
lextend = 1000
cut_off = 1
input = None
output = None
fullpath = os.path.realpath(__file__)
filename = sys.argv[0].split('/')[-1]
path = fullpath.split(filename)[0]
annotation = path+'hg19.ucsc.gene.txt'
try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:n:l:r:g:h')
except getopt.GetoptError as err:
	print str(err)
	usage()
	sys.exit(2)
for o,a in opts:
	if o == '-n':
		cut_off = int(a)
	elif o == '-l':
		lextend = int(a)
	elif o == '-i':
		input = a
	elif o == '-r':
		lread = int(a)
	elif o == '-o':
		output = a
	elif o == '-g':
		annotaton = a
	elif o == '-h':
		usage()
		sys.exit()
	else:
		assert False, "unhandled option"
if not input or not output:
	usage()
	sys.exit(2)

pair_set = set() 
cluster = []
ifp = open(input)
ofp = open(input+'.sv.tmp.txt','w')
for i,line in enumerate(ifp):
	item = line.rstrip().split()
	if item[0] in pair_set:
		continue
	else:
		pair_set.add(item[0])

	#get cluster for pairs
	iv1 = HTSeq.GenomicInterval(item[2],int(item[3])-lextend,int(item[3])+lread+lextend,'.')
	if 'chr' in item[6]:
		iv2 = HTSeq.GenomicInterval(item[6],int(item[7])-lextend,int(item[7])+lread+lextend,'.')
	elif item[6] == '=':
		iv2 = HTSeq.GenomicInterval(item[2],int(item[7])-lextend,int(item[7])+lread+lextend,'.')
	else:
		iv2 = HTSeq.GenomicInterval(item[2],int(item[6])-lextend,int(item[6])+lread+lextend,'.')
	merged = False
	for each in cluster:
		if each[0].overlaps(iv1):
			if each[1].overlaps(iv2):
				extend_interval(each[0],iv1)
				extend_interval(each[1],iv2)
				each[2] += 1
				merged = True
				break
	if not merged:
		cluster.append([iv1,iv2,1])

for each in cluster:
	if each[2] >= cut_off:
		print >> ofp, input+'\t'+str(each[2])+'\t'+each[0].chrom+'\t'+str(each[0].start+lextend)+'\t'+str(each[0].end-lextend)+'\t'+each[1].chrom+'\t'+str(each[1].start+lextend)+'\t'+str(each[1].end-lextend)
ifp.close()
ofp.close()

# adding genomic annotation information to detected sv file.
ifp = open(input+'.sv.tmp.txt')
ofp = open(output,'a')
print >> ofp, 'type\tdiscordant_pairs\tleft_chr\tleft_start\tleft_end\tright_chr\tright_start\tright_end\tleft_ann\tright_ann\tgene_fusion'
ofp1 = open(input+'.left.bed','w')
ofp2 = open(input+'.right.bed','w')
for line in ifp:
	items = line.rstrip().split()
	p = re.match('chr(.*)',items[2])
	if p.group(1) == 'X':
		items[2] = '23'
	elif p.group(1) == 'Y':
		items[2] = '24'
	else:
		items[2] = p.group(1)
	p = re.match('chr(.*)',items[5])
	if p.group(1) == 'X':
		items[5] = '23'
	elif p.group(1) == 'Y':
		items[5] = '24'
	else:
		items[5] = p.group(1)
	print >> ofp1,'\t'.join(items[2:5])
	print >> ofp2,'\t'.join(items[5:8])
ifp.close()
ofp1.close()
ofp2.close()
os.system('perl '+path+'addcatann.pl human 10 2 2 '+input+'.left.bed 0 '+annotation+' '+input+'.left.out')
os.system('perl '+path+'addcatann.pl human 10 2 2 '+input+'.right.bed 0 '+annotation+' '+input+'.right.out')
ifp = open(input+'.sv.tmp.txt')
ifp1 = open(input+'.left.out')
ifp2 = open(input+'.right.out')
while True:
	line = ifp.readline().rstrip()
	if not line:
		break
	line1 = ifp1.readline().rstrip().split('\t')
	line2 = ifp2.readline().rstrip().split('\t')
	l_ann = ','.join(line1[3:])
	r_ann = ','.join(line2[3:])
	if 'intergenic' not in l_ann and 'intergenic' not in r_ann:
		fusion = 'YES'
	else:
		fusion = 'NO'
	print >> ofp, line+'\t'+l_ann+'\t'+r_ann+'\t'+fusion
ifp.close()
ifp1.close()
ifp2.close()
ofp.close()
os.system('rm '+input+'.sv.tmp.txt')
os.system('rm '+input+'.left*')
os.system('rm '+input+'.right*')

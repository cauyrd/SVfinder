# usage: python classify_reads.py samfile_with_dup samfile_no_dup max_insert_size
#history: Wed Feb 22 13:19:27 EST 2012
import sys
import numpy as np
from scipy import stats
def get_uniq_pairs_list(ifp):
	unique_discordant_pairs = set()
	while True:
		line1 = ifp.readline().rstrip()
		if not line1:
			break
		if line1[0] == '@':
			continue
		line2 = ifp.readline().rstrip()
		items = line1.split()
		flag_num = int(items[1])
		if not (flag_num&0x4) and not(flag_num&0x8):
			if 'XT:A:U' in line1 and 'XT:A:U' in line2 and not flag_num&0x2:
				unique_discordant_pairs.add(items[0])
	return unique_discordant_pairs

if __name__ == '__main__':
	ifp = open(sys.argv[1])
	ofp_difchr = open('diffchr_pairs.txt','w')
	ofp_FR = open('discordantFR.txt','w')
	ofp_RF = open('discordantRF.txt','w')
	ofp_FF = open('discordantFF.txt','w')
	ofp_RR = open('discordantRR.txt','w')
	ofp_longisize = open('longisize.txt','w')
	ofp_uniq_left_re = open('uniq_left_read_remains.txt','w')
	ofp_uniq_right_re = open('uniq_right_read_remains.txt','w')
	ofp_proper = open('uniq_propely_mapped_reads.txt','w')
	unique_discordant_pairs = get_uniq_pairs_list(ifp)
	ifp.close()
	ifp = open(sys.argv[2])
	summary = {} 
	isize_set = []
	for line in ifp:
		if line[0] == '@':
			continue
		items = line.rstrip().split()
		flag_num = int(items[1])
		# select all mapped pairs
		if  not (flag_num & 0x4) and not (flag_num & 0x8):
			try:
				summary['mapped_pairs'] += 1 # count the number of reads (if count the number of pairs, divide this number by 2)
			except KeyError:
				summary['mapped_pairs'] = 1 
			# select proper pairs
			if flag_num & 0x2:
				print >> ofp_proper, line.rstrip()
				try:
					summary['proper_pairs'] += 1 
				except KeyError: 
					summary['proper_pairs'] = 1 
				if int(items[8])>0:
					isize_set.append(int(items[8]))
			# select discordant pairs with both ends uniquely mapped
			elif items[0] in unique_discordant_pairs: 
				if items[6] != '=': 
					print >> ofp_difchr, line.rstrip()
					try:
						summary['diffchr_pairs'] += 1 
					except KeyError:
						summary['diffchr_pairs'] = 1 
						
				elif (abs(int(items[8])) > int(sys.argv[3])):
					print >> ofp_longisize, line.rstrip()
					try:
						summary['longisize_pairs'] += 1 
					except KeyError:
						summary['longisize_pairs'] = 1 
				elif (flag_num & 0x10) and (flag_num & 0x20):
					print >> ofp_RR, line.rstrip()
					try:
						summary['discordantRR_pairs'] += 1 
					except KeyError: 
						summary['discordantRR_pairs'] = 1 
				elif not (flag_num & 0x10) and not (flag_num & 0x20):
					print >> ofp_FF, line.rstrip()
					try:
						summary['discordantFF_pairs'] += 1 
					except KeyError:
						summary['discordantFF_pairs'] = 1
				elif ((flag_num & 0x10) and (flag_num & 0x40)) or (not(flag_num & 0x10) and (flag_num & 0x80)):
					print >> ofp_RF, line.rstrip()
					try:
						summary['discordantRF_pairs'] += 1
					except KeyError:
						summary['discordantRF_pairs'] = 1 
				else:
					print >> ofp_FR, line.rstrip()
					try:
						summary['discordantFR_minus_pairs'] += 1 
					except KeyError:
						summary['discordantFR_minus_pairs'] = 1 
			else:
				try:
					summary['repetitive_pairs'] += 1
				except KeyError:
					summary['repetitive_pairs'] = 1
				if flag_num&0x40 and 'XT:A:U' in line:
					print >> ofp_uniq_left_re, line.rstrip()
					try:
						summary['uniquely_mapped_left_reads_remains'] += 1
					except KeyError:
						summary['uniquely_mapped_left_reads_remains'] = 1
				elif flag_num&0x80 and 'XT:A:U' in line:
					print >> ofp_uniq_right_re, line.rstrip()
					try:
						summary['uniquely_mapped_right_reads_remains'] += 1
					except KeyError:
						summary['uniquely_mapped_right_reads_remains'] = 1
		elif (flag_num & 0x8) and not (flag_num & 0x4):
			try:
				summary['singletons'] += 1 
			except KeyError:
				summary['singletons'] = 1
			if flag_num&0x40 and 'XT:A:U' in line:
				print >> ofp_uniq_left_re, line.rstrip()
				try:
					summary['uniquely_mapped_left_reads_remains'] += 1
				except KeyError:
					summary['uniquely_mapped_left_reads_remains'] = 1
			elif flag_num&0x80 and 'XT:A:U' in line:
				print >> ofp_uniq_right_re, line.rstrip()
				try:
					summary['uniquely_mapped_right_reads_remains'] += 1
				except KeyError:
					summary['uniquely_mapped_right_reads_remains'] = 1
		if not (flag_num & 0x4) and flag_num&0x40 and 'XT:A:U' in line:
			try:
				summary['uniquely_mapped_left_reads'] += 1
			except KeyError:
				summary['uniquely_mapped_left_reads'] = 1
		elif not (flag_num & 0x4) and flag_num&0x80 and 'XT:A:U' in line:
			try:
				summary['uniquely_mapped_right_reads'] += 1
			except KeyError:
				summary['uniquely_mapped_right_reads'] = 1
	isize_set = np.array(isize_set)
	bottom = stats.scoreatpercentile(isize_set,5)
	top = stats.scoreatpercentile(isize_set,95)
	isize_set = isize_set[(isize_set>bottom)*(isize_set<top)]
	mu = np.mean(isize_set)
	std = np.std(isize_set)
	for mykey in sorted(summary):
		print mykey+'\t'+str(summary[mykey])
	print 'mean:',mu
	print 'std:',std

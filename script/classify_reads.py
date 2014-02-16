# usage: python classify_reads.py mapped_read.sam
# history: 
#	1.Wed Feb 22 13:19:27 EST 2012
#	2.Wed Jan 15 16:30:15 EST 2014
import sys
import numpy as np
from scipy import stats
def isize_estimation(ifp):
	isize_set = []
	for line in ifp:
		if line[0] == '@':
			continue
		items = line.rstrip().split()
		flag_num = int(items[1])
		if not (flag_num&0x4) and not (flag_num&0x8):
			try:
				isize = int(items[8])
			except:
				isize = 0
			if flag_num&0x2 and isize > 0:
				isize_set.append(isize)
	isize_set = np.array(isize_set)
	bottom = stats.scoreatpercentile(isize_set,5)
	top = stats.scoreatpercentile(isize_set,95)
	isize_set = isize_set[(isize_set>bottom)*(isize_set<top)]
	mu = int(np.mean(isize_set))
	std = int(np.std(isize_set))
	max_isize = mu+3*std
	min_isize = np.abs(mu-3*std)
	print 'mean iszie:',mu,'std isize:',std,'max isize:',max_isize,'min isize:',min_isize
	return max_isize, min_isize

if __name__ == '__main__':
	ifp = open(sys.argv[1])
	ofp_difchr = open('diffchr_pairs.svReads.sam','w')
	ofp_FR = open('discordantFR.svReads.sam','w')
	ofp_RF = open('discordantRF.svReads.sam','w')
	ofp_FF = open('discordantFF.svReads.sam','w')
	ofp_RR = open('discordantRR.svReads.sam','w')
	ofp_longisize = open('longisize.svReads.sam','w')
	ofp_shortisize = open('shortisize.svReads.sam','w')
	ofp_singleton = open('singleton.svReads.sam','w')
	ofp_unmapped = open('unmapped.svReads.sam','w')
	max_isize, min_isize = isize_estimation(ifp)
	ifp.close()
	ifp = open(sys.argv[1])
	summary = {} 
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
				try:
					summary['proper_pairs'] += 1 
				except KeyError: 
					summary['proper_pairs'] = 1 
			# select discordant pairs with both ends uniquely mapped
			elif items[6] != '=': 
				print >> ofp_difchr, line.rstrip()
				try:
					summary['diffchr_pairs'] += 1 
				except KeyError:
					summary['diffchr_pairs'] = 1 
						
			elif abs(int(items[8])) > max_isize:
				print >> ofp_longisize, line.rstrip()
				try:
					summary['longisize_pairs'] += 1 
				except KeyError:
					summary['longisize_pairs'] = 1 
			elif abs(int(items[8])) <  min_isize:
				print >> ofp_shortisize, line.rstrip()
				try:
					summary['shortisize_pairs'] += 1 
				except KeyError:
					summary['shortisize_pairs'] = 1 
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
		elif (flag_num & 0x8) and not (flag_num & 0x4):
			print >> ofp_singleton, line.rstrip()
			try:
				summary['singletons'] += 1 
			except KeyError:
				summary['singletons'] = 1
		else:
			print >> ofp_unmapped, line.rstrip()
			try:
				summary['unmapped'] += 1 
			except KeyError:
				summary['unmapped'] = 1
	for mykey in sorted(summary):
		print mykey+'\t'+str(summary[mykey])

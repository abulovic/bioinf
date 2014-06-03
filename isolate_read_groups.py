#!/usr/bin/python
###############################################################
# This script is intented for isolating different read groups #
# from the bwa mem alignment. 								  #
# The assumption is that if there are multiple HIV strans     #
# present in the sequencing sample, there will be multiple    #
# read goups which will differ in certain places and will all #
# be mutually alike.                                          #
###############################################################

import sys
from collections import defaultdict

import pysam


def process_read(read, pos_distr, interesting_len=200):
	start = read.aend - read.alen
	end = read.aend
	for i in xrange(start, end - interesting_len):
		subseq = read.query[i - start: i - start + interesting_len]
		if subseq not in pos_distr[i]:
			pos_distr[i][subseq] = 1
		else:
			pos_distr[i][subseq] += 1


def main():
	if len(sys.argv) < 3:
		print '\nUsage: python isolate_read_groups.py [INPUT_BAM_FILE] [OUTPUT_DIR]\n'
		sys.exit(-1)
	in_fname = sys.argv[1]
	output_dir = sys.argv[2]

	# dictionary for mapping the positions of the alignment
	# to the read motifs aligned to that region of the ref genome
	pos_distribution = defaultdict(dict)

	# Start filling the dict by processing read by read
	# all the positions where the reads maps to
	print 'Loading reads coverage...'
	bam_in = pysam.Samfile(in_fname, 'rb')
	for read in bam_in.fetch():
		process_read(read, pos_distribution)
	print 'Reads coverage laoded.'

	# Isolate the positions where there are two or more
	# groups of reads that map equaly well to some 
	# position on the reference genome
	ambiguous = 0
	for pos, covs_dict in pos_distribution.iteritems():
		covs = sorted(covs_dict.values(), reverse=True)
		if len(covs) > 1 and float(covs[1]) / covs[0] > 0.2:
			ambiguous += 1
			print pos, float(covs[1]) / covs[0]


if __name__ == '__main__':
	main()
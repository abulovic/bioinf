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

import pylab as plab
from numpy import NaN, Inf, arange, isscalar, asarray, array
from scipy.signal import argrelextrema
import pysam

def peakdet(v, delta, x = None):
	maxtab = []
	mintab = []
	if x is None:
		x = arange(len(v))
	v = asarray(v)
	if len(v) != len(x):
		sys.exit('Input vectors v and x must have same length')
	if not isscalar(delta):
		sys.exit('Input argument delta must be a scalar')
	if delta <= 0:
		sys.exit('Input argument delta must be positive')
	mn, mx = Inf, -Inf
	mnpos, mxpos = NaN, NaN

	lookformax = True
	for i in arange(len(v)):
		this = v[i]
		if this > mx:
			mx = this
			mxpos = x[i]
		if this < mn:
			mn = this
			mnpos = x[i]
		if lookformax:
			if this < mx-delta:
				maxtab.append((mxpos, mx))
				mn = this
				mnpos = x[i]
				lookformax = False
		else:
			if this > mn+delta:
				mintab.append((mnpos, mn))
				mx = this
				mxpos = x[i]
				lookformax = True

	return array(maxtab), array(mintab)



def process_read(read, pos_distr, interesting_len=100):
	start = read.aend - read.alen
	end = read.aend
	for i in xrange(start, end - interesting_len):
		subseq = read.query[i - start: i - start + interesting_len]
		if subseq not in pos_distr[i]:
			pos_distr[i][subseq] = 1
		else:
			pos_distr[i][subseq] += 1

def plot_multimodal_read_distr(multimodal):
	x = range(0, max(multimodal.keys()) + 1)
	y = array([0.] * len(x))
	for key, val in multimodal.items():
		y[key] = val
	maxima, gl = peakdet(y, 0.3)
	lmx = []
	lmy = []
	for s in maxima:
		lmx.append(int(s[0]))
		lmy.append(s[1])
	plab.plot(x, y)
	plab.hold(True)
	plab.plot(lmx, lmy, 'ro')
	plab.axis([0, len(y), 0, 2])
	plab.show()


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
	multimodal = {}
	for pos, covs_dict in pos_distribution.iteritems():
		covs = sorted(covs_dict.values(), reverse=True)
		if len(covs) > 1 and float(covs[1]) / covs[0] > 0.2:
			ambiguous += 1
			multimodal[pos] = float(covs[1]) / covs[0]

	plot_multimodal_read_distr(multimodal)

if __name__ == '__main__':
	main()
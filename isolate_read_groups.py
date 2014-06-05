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
import math
from collections import defaultdict

import pylab as plab
from numpy import NaN, Inf, arange, isscalar, asarray, array
from scipy.signal import argrelextrema
import pysam

WINDOW_SIZE = 100
MIN_ALNS = 50

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



def process_read(read, pos_distr, interesting_len=WINDOW_SIZE):
	start = read.aend - read.alen
	end = read.aend
	for i in xrange(start, end - interesting_len):
		subseq = read.query[i - start: i - start + interesting_len]
		if subseq not in pos_distr[i]:
			pos_distr[i][subseq] = 1
		else:
			pos_distr[i][subseq] += 1

def get_signal(multimodal):
	x = range(0, max(multimodal.keys()) + 1)
	y = array([0.] * len(x))
	for key, val in multimodal.items():
		y[key] = val
	return x, y

def get_local_maxima(sig):
	maxima, gl = peakdet(sig, 0.3)
	lmx = []
	lmy = []
	for s in maxima:
		lmx.append(int(s[0]))
		lmy.append(s[1])
	return lmx, lmy

def calc_entropy(probs):
	return - sum(map(lambda p: p * math.log(p, 2), probs))

def get_entropy(probs_distribution, scaled=False):
	entropy = [0] * (max(probs_distribution.keys()) + 1)
	for loc, distr in probs_distribution.iteritems():
		total_reads = sum(distr.values())
		probs = map(lambda n: float(n) / total_reads, distr.values())
		entropy[loc] = calc_entropy(probs)
		if scaled:
			entropy[loc] *= total_reads
	max_entropy = max(entropy)
	entropy = map(lambda e: e / max_entropy, entropy)
	return entropy


def plot_multimodal_read_distr(x, y, lmx, lmy):
	plab.figure()
	plab.title('Genome locations with 2-modal read distributions, window length = %d' % WINDOW_SIZE)
	plab.plot(x, y)
	plab.hold(True)
	plab.plot(lmx, lmy, 'ro')
	plab.axis([0, len(y), 0, 2])


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
		if len(covs) > 1 and float(covs[1]) / covs[0] > 0.2 and sum(covs) > MIN_ALNS:
			ambiguous += 1
			multimodal[pos] = float(covs[1]) / covs[0]

	x, y = get_signal(multimodal)
	lmx, lmy = get_local_maxima(y)
	plot_multimodal_read_distr(x, y, lmx, lmy)

	# We calculate the scaled entropy. 
	# The higher the entropy, the more versatile the reads aligned
	# to that location. 
	# We scale the entropy with the number of alignments to that position,
	# since we also want it to indicate whether there were many or few
	# alignments involved.
	entropy = get_entropy(pos_distribution, scaled=True)
	plab.plot(entropy, 'g')

	for loc in lmx:
		print "LOC %4d: ENTROPY: %1.4f, TOTAL READS: %d, SEC/FIRST: %1.4f" % (loc, entropy[loc], sum(pos_distribution[loc].values()), multimodal[loc])

	plab.show()

if __name__ == '__main__':
	main()
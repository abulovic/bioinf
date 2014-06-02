#!/usr/bin/env python

import sys

import numpy as np
import pylab as P
import pysam as ps


bam_fname = sys.argv[1]
sam_file = ps.Samfile(bam_fname)

mapq = map(lambda read: read.mapq, iter(sam_file.fetch()))
rlen = map(lambda read: read.rlen, iter(sam_file.fetch()))

P.figure()
n, bins, patches = P.hist(mapq, 60)
P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

P.figure()
nbins = max(rlen) - min(rlen) + 1
n, bins, patches = P.hist(rlen, nbins)
P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)


P.show()
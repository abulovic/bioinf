import sys
from contextlib import nested

import pysam

in_fname = sys.argv[1]

in_file = pysam.Samfile(in_fname, 'rb')
for read in in_file.fetch():
	if read.aend > 7000:
		import pdb
		pdb.set_trace()
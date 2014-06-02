import sys
from contextlib import nested

import pysam

in_fname = sys.argv[1]
out_fname = sys.argv[2]

in_file = pysam.Samfile(in_fname, 'rb')
out_file = pysam.Samfile(out_fname, 'wb', template=in_file)
num_reads = 0

for read in in_file.fetch():
	if read.mapq < 60L:
		num_reads += 1
		out_file.write(read)

print 'Number of imperfect reads: ', num_reads
in_file.close()
out_file.close()
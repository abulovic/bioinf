import sys

import pysam


perfect = set()
imperfect = set()

bam_fname = sys.argv[1]
fin = pysam.Samfile(bam_fname, 'rb')

for read in fin.fetch():
	if read.mapq == 60L:
		perfect.add(read.qname)
	else:
		imperfect.add(read.qname)
fin.close()

unused_fastq = sys.argv[2]
unused = 0
with open(unused_fastq, 'r') as fin:
	while True:
		l = fin.readline().strip()
		if l == '':
			break
		if l.startswith('@'):
			unused += 1

print 'Perfect:', len(perfect)
print 'Imperfect', len(imperfect)
print 'Unused', unused
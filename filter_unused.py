import sys
from contextlib import nested

import pysam

def load_id(f):
	rid = f.readline().strip()
	if rid == '':
		return ''
	for i in range(0,3):
		f.readline()
	return rid.split()[0][1:]

def load_entry(f):
	rid = f.readline().strip()
	if rid == '':
		return '', '', ''
	seq = f.readline().strip()
	f.readline()
	phred = f.readline().strip()
	return rid, seq, phred

in1 = sys.argv[1]
in2 = sys.argv[2]
bam = sys.argv[3]
out = sys.argv[4]
out1 = '%s1.fastq' % out
out2 = '%s2.fastq' % out

fastq1_ids = set()
fastq2_ids = set()
aln_ids = set()

bam_file = pysam.Samfile(bam, 'rb')
for read in bam_file.fetch():
	aln_ids.add(read.qname)

with nested(open(in1), open(in2)) as (f1, f2):
	while True:
		rid1 = load_id(f1)
		if rid1 == '':
			break
		rid2 = load_id(f2)
		fastq1_ids.add(rid1)
		fastq1_ids.add(rid2)

nonaligned = (fastq1_ids | fastq2_ids) - aln_ids

with nested(open(in1, 'r'), open(out1, 'w')) as (fin, fout):
	while True:
		header, seq, phred = load_entry(fin)
		if header == '':
			break
		rid = header.split()[0][1:]
		if rid in nonaligned:
			fout.write('%s\n%s\n+\n%s\n' % (header, seq, phred))

with nested(open(in2, 'r'), open(out2, 'w')) as (fin, fout):
	while True:
		header, seq, phred = load_entry(fin)
		if header == '':
			break
		rid = header.split()[0][1:]
		if rid in nonaligned:
			fout.write('%s\n%s\n+\n%s\n' % (header, seq, phred))


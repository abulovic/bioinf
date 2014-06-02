import sys
from contextlib import nested

def main():

	seq_data1 = dict()
	seq_data2 = dict()

	fname1 = sys.argv[1]
	fname2 = sys.argv[2]

	with open(fname1, 'r') as f1:
		while True:
			seq_id = f1.readline().strip()
			if seq_id == '':
				break
			seq = f1.readline().strip()
			f1.readline()
			quality = f1.readline().strip()

			seq_unique_id = seq_id.split()[0]
			seq_data1[seq_unique_id] = (seq_id, seq, quality)

	with open(fname2, 'r') as f2:
		while True:
			seq_id = f2.readline().strip()
			if seq_id == '':
				break
			seq = f2.readline().strip()
			f2.readline()
			quality = f2.readline().strip()

			seq_unique_id = seq_id.split()[0]
			seq_data2[seq_unique_id] = (seq_id, seq, quality)

	seqs1 = set(seq_data1.keys())
	seqs2 = set(seq_data2.keys())

	seqs_intersection = seqs1 & seqs2
	out_fname1 = '%s.clean.fastq' % fname1.split('.')[0]
	out_fname2 = '%s.clean.fastq' % fname2.split('.')[0]
	with nested(open(out_fname1, 'w'), open(out_fname2, 'w')) as (f1, f2):
		for seq in seqs_intersection:
			seq_id1, seq1, quality1 = seq_data1[seq]
			f1.write('%s\n%s\n%s\n%s\n' % (seq_id1, seq1, '+', quality1))
			seq_id2, seq2, quality2 = seq_data2[seq]
			f2.write('%s\n%s\n%s\n%s\n' % (seq_id2, seq2, '+', quality2))

if __name__ == '__main__':
	main()
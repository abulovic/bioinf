import sys
from contextlib import nested

def main():
	input_fname = sys.argv[1]
	output_fname = sys.argv[2]

	seq_set = set()
	total_seqs = 0
	double_seqs = 0
	seqs_with_n = 0

	with nested(open(input_fname, 'r'), open(output_fname, 'w')) as (in_file, out_file):
		while True:
			seq_id = in_file.readline().strip()
			if seq_id == '':
				break
			seq = in_file.readline().strip()
			in_file.readline()
			seq_quality = in_file.readline().strip()

			total_seqs += 1

			if seq in seq_set:
				double_seqs += 1
				continue
			if 'N' in seq:
				seqs_with_n += 1
				continue

			seq_set.add(seq)
			out_file.write('%s\n%s\n%s\n%s\n' % (seq_id, seq, '+', seq_quality))

		print 'Number of seqs in the original file: ', total_seqs
		print 'Number of double sequences: ', double_seqs
		print 'Number of seqs with Ns:     ', seqs_with_n


if __name__ == '__main__':
	main()
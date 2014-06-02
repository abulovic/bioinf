# python ~/Software/src/HIV/scripts/analyze_similarity.py NIBSC_stats.txt `ls NI*.aln`
import sys

def load_seqs(f):
	lines = f.readlines()
	for idx, line in enumerate(lines[1:]):
		if line.startswith('>'):
			break
	idx += 1
	first_seq = ''.join(map(lambda s: s.strip(), lines[1:idx]))
	second_seq = ''.join(map(lambda s: s.strip(), lines[idx+1:]))
	return first_seq, second_seq

stats_fname = sys.argv[1]
with open(stats_fname, 'w') as stats_file:
	for fname in sys.argv[2:]:
		print fname
		with open(fname) as f:
			ref_genome, assembly = load_seqs(f)
			total = len(ref_genome)
			insertions = ref_genome.count('-')
			deletions = assembly.count('-')
			mismatches = 0
			for a, b in zip(ref_genome, assembly):
				if a != b and a != '-':
					mismatches += 1
			identical = total - insertions - deletions - mismatches
		stats_file.write('%s %f %f %f %f\n' % (fname.split('.')[0],
											 float(identical) / (total - insertions - deletions) * 100.,
											 float(mismatches) / (total - insertions - deletions) * 100.,
											 float(deletions) / total * 100.,
											 float(insertions) / total * 100.))

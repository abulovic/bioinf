import sys
from contextlib import nested

fqname = sys.argv[1]
faname = '.'.join(fqname.split('.')[:-1]) + '.fa'

cnt = -1
with nested(open(fqname), open(faname, 'w')) as (fin, fout):
	for l in fin:
		cnt += 1
		l = l.strip()
		if cnt % 4 == 0:
			l = '>' + l[1:]
		if cnt % 4 in (2, 3):
			continue
		fout.write('%s\n' % l)

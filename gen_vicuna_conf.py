#!/usr/bin/python
import sys

if len(sys.argv) < 4:
	print "Usage: python gen_vicuna_conf.py "\
	"[INPUT_DIR] [OUTPUT_DIR] [VICUNA_CONF_DIR]"
	sys.exit(-1)

input_dir = sys.argv[1]
output_dir = sys.argv[2]
vicuna_conf_dir = sys.argv[3]

vconf_fname = '%s/vicuna.conf' % vicuna_conf_dir
with open(vconf_fname, 'w') as vconf_file:
	vconf_file.write('maxOverhangSize	2\n')
	vconf_file.write('Divergence	8\n')
	vconf_file.write('max_read_overhang	8\n')
	vconf_file.write('max_contig_overhang	10\n')
	vconf_file.write('pFqDir	%s/\n' % input_dir)
	vconf_file.write('batchSize	1000000\n')
	vconf_file.write('LibSizeLowerBound	100\n')
	vconf_file.write('LibSizeUpperBound	800\n')
	vconf_file.write('min_output_contig_len	200\n')
	vconf_file.write('outputDIR	%s/\n' % output_dir)

vanalasys_conf_fname = '%s/vanalasys.conf' % vicuna_conf_dir
with open(vanalasys_conf_fname, 'w') as vanconf_file:
	vanconf_file.write('pFqDir	%s/\n' % input_dir)
	vanconf_file.write('trim_log_file %s/trim.log\n' % output_dir)
	vanconf_file.write('aln_file %s/contig.align\n' % output_dir)
	vanconf_file.write('outputDIR %s/\n' % output_dir)
	vanconf_file.write('num_region 1\n')
	vanconf_file.write('2	0	500\n')

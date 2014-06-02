#!/bin/bash

R1_FILE=$1
R2_FILE=$2
OUTPUT_DIR=$3

# Prepare reads for assembly
# Let's first inspect the quality of data (PREQC)
# pe-mode 1: if read is discarded, its pair is discarded also
sga preprocess --verbose --pe-mode 1 $R1_FILE $R1_FILE > $OUTPUT_DIR/sga.reads.fastq

# Index the reads in READSFILE using a suffixarray/bwt
sga index --verbose -a ropebwt -t 8 $OUTPUT_DIR/sga.reads.fastq 

# Perform pre-assembly quality checks
# sga preqc -t 8 $OUTPUT_DIR/sga.reads.fastq > $OUTPUT_DIR/sga.reads.preqc
# !!! For some reason, sga preqc gets stuck at "Loading FM-index of sga.reads.fastq"
# !!! Update, OK not stuck, just impossibly long... :D
# sga-preqc-report.py $OUTPUT_DIR/sga.reads.preqc

# Correct sequencing errors in all the reads
sga correct --verbose --learn -o $OUTPUT_DIR/sga.reads.corr.fastq $OUTPUT_DIR/sga.reads.fastq

# Index the corrected data
sga index --verbose -a ropebwt -t 8 $OUTPUT_DIR/sga.reads.corr.fastq

# Remove the exact match duplicates and reads with low frequency k-mers
sga filter -x 2 -t 8 $OUTPUT_DIR/sga.reads.corr.fastq

# Compute the structure of the string graph using 
# the corrected and filtered reads
sga overlap -t 8 $OUTPUT_DIR/sga.reads.corr.filter.pass.fa
# sga overlap -t 8 $OUTPUT_DIR/sga.reads.fastq

# Create contigs from the Assembly Graph ASQG file
sga assemble --verbose -m 30 --min-branch-length 400 -o primary $OUTPUT_DIR/sga.reads.corr.filter.pass.asqg.gz 
# sga assemble --verbose -m 30 --min-branch-length 400 -o primary $OUTPUT_DIR/sga.reads.asqg.gz 

# Align the reads to contigs
bwa index primary-contigs.fa
bwa aln -t 8  primary-contigs.fa $R1_FILE > $R1_FILE.sai
bwa aln -t 8  primary-contigs.fa $R2_FILE > $R2_FILE.sai
bwa sampe primary-contigs.fa $R1_FILE.sai $R2_FILE.sai $R1_FILE $R2_FILE | samtools view -Sb - > libPE.bam

# Convert the BAM file into a set of contig-contig distance estimates
sga-bam2de.pl -n 10 -m 250 --prefix libPE libPE.bam

# Compute copy number estimates of the contigs
sga-astat.py -m 250 libPE.bam > libPE.astat

# Build the scaffolds
sga scaffold -m 250 -a libPE.astat -o scaffolds.scaf --pe libPE.de primary-contigs.fa


#!/bin/bash
#####################################################################
################ RUN BWA ON PAIRED READS ############################

while getopts 1:2:r:o: opt; do
  # echo $opt
  case $opt in
  1)
      IN1=$OPTARG
      ;;
  2)
      IN2=$OPTARG
      ;;
  r)
      REF_GENOME=$OPTARG
      ;;
  o)
	  OUT=$OPTARG
	  ;;
  esac
done

if [[ -z $IN1 || -z $IN2 || -z $REF_GENOME || -z $OUT ]]; then
	echo 'Not all required arguments supplied.'
	echo 'Run the program with following parameters:'
	echo -e '\t-1: first fastq reads file'
	echo -e '\t-2: second fastq reads file'
	echo -e '\t-r: reference genome fasta'
	echo -e '\t-o: output file prefix (no file extension)'
	exit
fi

OUT_SAM=$OUT.sam
OUT_BAM=$OUT.bam
OUT_SORTED=$OUT.sorted
OUT_SORTED_BAM=$OUT.sorted.bam
OUT_SORTED_BAI=$OUT.sorted.bai
OUT_IMPERFECT_BAM=$OUT.imperfect.bam
OUT_IMPERFECT_BAI=$OUT.imperfect.bai

# Path to the bwa executable (in case of multiple versions installed)
BWA=bwa2

# Check if the right version of the bwa software is installed
BWA_VERSION=`$BWA 2>/dev/stdout | grep 'Version' | sed 's/Version: 0.\([0-9]\).*/\1/'`
if [ $BWA_VERSION -lt 7 ]; then
	echo 'BWA version 0.7.x requred for usage of bwa-mem algorithm.'
	exit
fi

# Indexing of the reference genome
echo 'Running bwa index $REF_GENOME...'
$BWA index $REF_GENOME

# 1. Run the bwa-mem
$BWA mem $REF_GENOME $IN1 $IN2 > $OUT_SAM
# 2. Convert .sam to .bam
samtools view -bS $OUT_SAM > $OUT_BAM
# 3. Sort the .bam file
samtools sort $OUT_BAM $OUT_SORTED
# 4. Index the .bam file
samtools index $OUT_SORTED_BAM $OUT_SORTED_BAI

# Post-alignment steps:
# 1. Create fastq with sequences which weren't used in the alignment
python filter_unused.py $IN1 $IN2 $OUT_SORTED_BAM $OUT.unused
# 2. Create a .bam file with seqs that have imperfect mapq score
python filter_imperfect_reads.py $OUT_SORTED_BAM $OUT_IMPERFECT_BAM
samtools index $OUT_IMPERFECT_BAM $OUT_IMPERFECT_BAI

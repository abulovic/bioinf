#!/bin/bash
#####################################################################
##### RUN PAIRWISE ALNS FOR DE NOVO GENOME ASSEMBLIES ###############

while getopts r:i:o:n: opt; do
  case $opt in
  r)
      REF_DIR=$OPTARG
      ;;
  i)
      INPUT_FILE=$OPTARG
      ;;
  o)
      OUTPUT_DIR=$OPTARG
      ;;
  n)
	  NAME=$OPTARG
	  ;;
  esac
done

if [[ -z $REF_DIR || -z $INPUT_FILE || -z $OUTPUT_DIR || -z $NAME ]]; then
	echo 'Not all required arguments supplied.'
	echo 'Run the program with following parameters:'
	echo -e '\t-i: Input file containing de novo assembly or contigs.'
	echo -e '\t-2: second fastq reads file'
	echo -e '\t-r: reference genome fasta'
	echo -e '\t-o: output file prefix (no file extension)'
	exit
fi

subtypes=(A A1 A2 B C D F1 F2 G H J K N)

for subtype in ${subtypes[*]}
do
	cat $REF_DIR/subtype_${subtype}.fasta $INPUT_FILE > $OUTPUT_DIR/${NAME}_4x_${subtype}.fasta
	mafft $OUTPUT_DIR/${NAME}_4x_${subtype}.fasta > $OUTPUT_DIR/${NAME}_4x_${subtype}.aln
done
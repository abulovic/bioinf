#!/bin/bash
#####################################################################
##### RUN MULTIPLE SEQ ALNS FOR DE NOVO GENOME ASSEMBLIES ###########

while getopts r:i:o:n: opt; do
  case $opt in
  r)
      REF_DIR=$OPTARG
      ;;
  i)
      INPUT_DIR=$OPTARG
      ;;
  o)
      OUTPUT_DIR=$OPTARG
      ;;
  n)
	  NAME=$OPTARG
	  ;;
  esac
done

if [[ -z $REF_DIR || -z $INPUT_DIR || -z $OUTPUT_DIR || -z $NAME ]]; then
	echo 'Not all required arguments supplied.'
	echo 'Run the program with following parameters:'
	echo -e '\t-i: Input (results) directory containing de novo assembly of idba, sga and vicuna.'
  echo -e '\t-r: reference genome fasta'
  echo -e '\t-o: output file prefix (no file extension)'
	echo -e '\t-n: new file name prefix'
	exit
fi

subtypes=(A A1 A2 B C D F1 F2 G H J K N)

for subtype in ${subtypes[*]}
do
	cat $REF_DIR/subtype_${subtype}.fasta $INPUT_DIR/idba/4x/all/consensus_fixed_assembly.fa $INPUT_DIR/sga/4x/all/consensus_fixed_assembly.fa $INPUT_DIR/vicuna/4x/consensus.fa> $OUTPUT_DIR/${NAME}_4x_${subtype}.fasta
	mafft $OUTPUT_DIR/${NAME}_4x_${subtype}.fasta > $OUTPUT_DIR/${NAME}_4x_${subtype}.aln

  cat $REF_DIR/subtype_${subtype}.fasta $INPUT_DIR/idba/F/all/consensus_fixed_assembly.fa $INPUT_DIR/sga/F/all/consensus_fixed_assembly.fa $INPUT_DIR/vicuna/F/consensus.fa> $OUTPUT_DIR/${NAME}_F_${subtype}.fasta
  mafft $OUTPUT_DIR/${NAME}_F_${subtype}.fasta > $OUTPUT_DIR/${NAME}_F_${subtype}.aln
done
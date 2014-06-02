#!/bin/bash

INPUT_REF=$1	# directory containing referent genomes
INPUT_DIR=$2	# directory containing assembly results
OUTPUT_DIR=$3	# directory where to write the files
NAME=$4

subtypes=(A A1 A2 B C D F1 F2 G H J K N)

for subtype in ${subtypes[*]}
do
	cat $INPUT_REF/subtype_${subtype}.fasta $INPUT_DIR/4x/consensus.fa > $OUTPUT_DIR/${NAME}_4x_${subtype}.fasta
	mafft $OUTPUT_DIR/${NAME}_4x_${subtype}.fasta > $OUTPUT_DIR/${NAME}_4x_${subtype}.aln
done
#!/bin/bash
#Script to take a bed file and fine the nearest gene.
#does this using bedtools- first sorts in place, then maps.
#Optimized for use on MARCC
#
set -e
BED_IN=$1
PROTEIN_CODING="/work-zfs/abattle4/ashton/reference_data/protein_coding_genes.bed"
OUT_FILE=$2
ml bedtools
echo "We first sort the bed file in place"
bedtools sort -i ${BED_IN} > ${BED_IN}.resorted.bed;
mv ${BED_IN}.resorted.bed ${BED_IN}

echo "We match each SNP to its closest gene"
bedtools closest -t first -a ${BED_IN} -b ${PROTEIN_CODING} -wb > ${OUT_FILE}

echo "Done!"
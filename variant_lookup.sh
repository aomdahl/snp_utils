#!/bin/bash
#Script to map RSIDs to chr locations

#RSIDs from 
##VCF downloaded from ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz on Feb 10, 2021
#Relevant columns extracted into 'rsid_chr_pos.txt', and header left in this file, by the following commands:

RSID=$1
OUT=$2
awk '(FNR == NR) {arr[$1];next} ($3 in arr) {print $1":"$2"\t"$3}' $RSID /work-zfs/abattle4/lab_data/hg19/variant_calls/rsid_chr_pos.txt > $OUT

#Now, we want them in the same order they started in.
awk '(FNR == NR) {arr[$2]=$1; next} ($1 in arr) {print $0}' $OUT $RSID > $OUT.tmp  && mv $OUT.tmp $OUT

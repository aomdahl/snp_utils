#!/bin/bash

ID=$1
OUT=$2
awk '(FNR == NR) {arr[$1];next} ($1":"$2 in arr) {print $1":"$2"\t"$3}' $ID /work-zfs/abattle4/lab_data/hg19/variant_calls/rsid_chr_pos.txt > $OUT
#Now get it back in the right order
awk '(FNR == NR) {arr[$1]=$2 ;next} ($1 in arr) {print $1"\t"arr[$1]}' $OUT $ID > ${OUT}.tmp
mv ${OUT}.tmp $OUT

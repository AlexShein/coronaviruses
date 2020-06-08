#!/bin/bash

# $1 - first arg - path to source fasta file with multiple sequences
# $2 - second arg - path to folder with processed fasta results
# example: src/split_fasta.sh all_cornaviruses_ncbi.fasta processed_fasta

RECORDS_COUNT=$(($(cat $1 | grep -c '>')-2))

csplit $1 '/^>.*/' "{$RECORDS_COUNT}"

for filename in $(ls xx*); do
    SEQ_NAME=$(
        head -n1 $filename | sed -e 's/^.*\.[0-9] \(.*\)$/\1/' -e 's/[ \/]/_/g' -e 's/[,>]//g'
        )
    echo "Moving $filename to $2/$SEQ_NAME.fasta"
    mv $filename $2/$SEQ_NAME.fasta
done;
echo "Done!"
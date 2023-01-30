#!/bin/bash
# conda activate biopiton
python getlengthfasta.py $1 $2
cut -f2 data/length_fasta.tsv | tail -n +2 | paste -d "\t" > data/metadata.tmp
paste -d "\t" data/JudiSeq-A25_mapping.txt data/metadata.tmp > data/metadata.tsv
rm -rf data/*.tmp data/length_*
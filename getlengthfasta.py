#!/usr/bin/env python

# USAGE: python getlengthfasta.py [fasta file] [outdir]
# fasta file: (required) fasta file in fasta format
# outdir: (optional, default: current dir) path to the dir where the output will be placed

import sys
from Bio import SeqIO
import pandas as pd

# Required functions
def parse_fasta(nombre):
    fasta_sequences = SeqIO.parse(nombre,'fasta')
    ids = []
    lengths = []
    for fasta in fasta_sequences:
        ids.append(fasta.id)
        lengths.append(len(fasta))
    df = pd.DataFrame(data={"seq_id": ids, "length": lengths})
    return df

# Arguments
raw_table = sys.argv[1]
if len(sys.argv) == 2:
    outdir = sys.argv[2]
else:
    outdir = "."

tabla = parse_fasta(raw_table)
with open(f"{outdir}/length_fasta.tsv", "w") as out:
    tabla.to_csv(out, sep = "\t", index = False, header = True)
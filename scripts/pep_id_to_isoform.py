#!/usr/bin/env python3

from Bio import SeqIO
import sys

def pep_to_isogroup(id):
    id_split = id.split("|")
    n = len(id_split)
    return "|".join(id_split[0:n-1])


if __name__ == "__main__":

    if len(sys.argv) != 1:
        sys.stderr.write("ERROR! Pipe me a fasta file with the transdecoder proteins (transdecoder.pep)!")

    if len(sys.argv) != 1:
        sys.stderr.write("ERROR! Pipe me a list of headers, one per line!")

    [sys.stdout.write(pep_to_isogroup(protein) + "\n") for protein in sys.stdin]
        

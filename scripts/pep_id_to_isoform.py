#!/usr/bin/env python3

from Bio import SeqIO
import sys

def pep_to_isogroup(id):
    return "|".join(id.split("|")[0:2])


if __name__ == "__main__":

    if len(sys.argv) != 1:
        sys.stderr.write("ERROR! Pipe me a fasta file with the transdecoder proteins (transdecoder.pep)!")

    if len(sys.argv) != 1:
        sys.stderr.write("ERROR! Pipe me a list of headers, one per line!")

    [sys.stdout.write(pep_to_isogroup(protein) + "\n") for protein in sys.stdin]
        

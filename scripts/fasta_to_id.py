#!/usr/bin/env python3

from Bio import SeqIO
import sys

if __name__ == "__main__":

    if len(sys.argv) != 1:
        sys.stderr.write("ERROR! Pipe me a fasta file!")

    [sys.stdout.write(record.id + "\n") 
        for record in SeqIO.parse(sys.stdin, "fasta")]
        

#!/usr/bin/env python3

# Read from stdin a series of ids of a trinity assembly and print the gene corresponding to it

import sys

def get_gene(header):
    """
    From a header of the form TRNNNN|gX_cY_iZ, get the TRNNNN|gX_cY and count how many times
    """
    header_splitted = header.split("_")
    return "_".join(header_splitted[0:2])


if __name__ == "__main__":

    if len(sys.argv) != 1:
        sys.stderr.write("ERROR! Pipe me a list of headers, one per line!")

    [sys.stdout.write(get_gene(isoform) + "\n") for isoform in sys.stdin]

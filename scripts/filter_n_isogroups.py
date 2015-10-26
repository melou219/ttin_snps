#!/usr/bin/env python3

from Bio import SeqIO
import sys

def read_headers(filename):
    """
    Read lines from a text file
    """
    headers = [line.strip() for line in open(filename, "r")]
    return headers



def get_gene(header):
    """
    From a header of the form TRNNNN|gX_cY_iZ, get the TRNNNN|gX_cY and count how many times.
    Works too on new Trinity headers with the form TRINITYNNNN_gX_cY_iZ
    """
    header_splitted = header.split("_")
    n = len(header_splitted)
    return "_".join(header_splitted[0:n-1])



def count_transcripts_per_gene(headers):
    """
    Number of times a Trinity gene is represented.
    """
    counts = {}
    for transcript in headers:
        gene = get_gene(transcript)
        if gene not in counts:
            counts[gene] = 1
        else:
            counts[gene] += 1
    return counts



def filter_multi_genes(headers):
    """
    Filter transcripts with genes appearing more than once in the header list.
    """
    counts = count_transcripts_per_gene(headers)
    one_isogroup = 0
    n_isogroup   = 0
    results = []
    for transcript in headers:
        if counts[get_gene(transcript)] == 1 :
            one_isogroup += 1
            results.append(transcript)
        else:
            n_isogroup += 1
    print("1-isogroups: %d\nN-isogroups: %d" % (one_isogroup, n_isogroup))
    return {transcript for transcript in headers if counts[get_gene(transcript)] == 1}



def write_headers(headers, filename_out):
    """
    Write headers to file_out, one per line.
    """
    with open(filename_out, "w") as f_out:
        [f_out.write(transcript + "\n") for transcript in headers]
    
    

if __name__ == "__main__":
    
    usage = "python3 filter_n_isogroups.py headers.in headers.out"
    
    if len(sys.argv) != 3:
        sys.exit("ERROR! Incorrect number of arguments \n" + usage)
    else:
        headers_in  = sys.argv[1]
        headers_out = sys.argv[2]
    
    headers_raw   = read_headers(headers_in)
    headers_final = filter_multi_genes(headers_raw)
    write_headers(headers_final, headers_out)

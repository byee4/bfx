#!/usr/bin/env python

"""
basic module to use as an example.
"""
import argparse
from Bio import SeqIO
import numpy as np

def n50(fasta):
    """
    Prints the n50 and other stats

    :param fasta:
    :return:
    """
    handle = open(fasta, "rU")
    all_lengths = []
    total_bases = 0
    total_contigs = 0
    for record in SeqIO.parse(handle, "fasta"):
        all_lengths.append(len(record.seq))
        total_bases += len(record.seq)
        total_contigs += 1
    all_lengths = sorted(all_lengths)
    print("Total bases: {}".format(total_bases))
    print("Total contigs: {}".format(total_contigs))
    print("Longest contig length: {}".format(max(all_lengths)))
    print("Shortest contig length: {}".format(min(all_lengths)))
    print("Mean contig length: {}".format(np.mean(all_lengths)))
    print("Median contig length: {}".format(np.median(all_lengths)))
    n = 0
    for l in all_lengths[::-1]:
        n += l
        if n > total_bases * (0.5):
            print("N50: {}".format(l))
            break

def main():
    """
    Main program.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fasta",
        required=True,
        default=None
    )
    args = parser.parse_args()
    fasta = args.fasta

    n50(fasta)

if __name__ == "__main__":
    main()

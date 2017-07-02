#!/usr/bin/env python

"""

Filters an input normalized peak BED6 file (log10 pvalue as the 'name' column,
log2 foldchange as the 'score' column).

"""
import argparse


def filter(input_file, output_file, l10p, l2fc):
    """
    Filters an input_file to have only peaks that pass l10p and l2fc cutoffs.

    :param input_file: basestring
    :param output_file: basestring
    :param l10p: float
    :param l2fc: float
    :return 0:
    """
    o = open(output_file, 'w')
    with open(input_file, 'r') as i:
        for line in i:
            line = line.split('\t')
            p = float(line[3])
            f = float(line[4])
            if p > l10p and f > l2fc:
                o.write('\t'.join(line))
    return 0

def main():
    """
    Main program.

    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        required=True,
        type=basestring,
    )
    parser.add_argument(
        "--output",
        required=False,
        type=basestring,
    )
    parser.add_argument(
        "--l10p",
        required=False,
        type=int,
        default=3
    )
    parser.add_argument(
        "--l2fc",
        required=False,
        type=int,
        default=3
    )
    args = parser.parse_args()

    l10p = args.l10p
    l2fc = args.l2fc
    i = args.input
    o = args.output if args.output is not None else i + '.filtered.bed'
    filter(i, o, l10p, l2fc)

if __name__ == "__main__":
    main()

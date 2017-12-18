import pandas as pd
import os
import sys
from argparse import ArgumentParser


def row_to_miso(row):
    if row['strand'] == '+':
        return '{}:{}:{}:{}@{}:{}:{}:{}@{}:{}:{}:{}\t{}'.format(
            row['chr'],int(row['upstreamES'])+1,row['upstreamEE'],row['strand'],
            row['chr'],int(row['exonStart_0base'])+1,row['exonEnd'],row['strand'],
            row['chr'],int(row['downstreamES'])+1,row['downstreamEE'],row['strand'],
            row['GeneID']
        )
    else:
        return '{}:{}:{}:{}@{}:{}:{}:{}@{}:{}:{}:{}\t{}'.format(
            row['chr'],int(row['downstreamES'])+1,row['downstreamEE'],row['strand'],
            row['chr'],int(row['exonStart_0base'])+1,row['exonEnd'],row['strand'],
            row['chr'],int(row['upstreamES'])+1,row['upstreamEE'],row['strand'],
            row['GeneID']
        )


def rmats_to_miso(input_file, output_file):
    """
    Writes a miso-like-formatted file from an rMATS JunctionCountsOnly one.

    Parameters
    ----------
    input_file : string
    output_file : string

    Returns
    -------
    None
    """

    df = pd.read_table(input_file)
    miso = df.apply(row_to_miso, axis=1)
    miso.to_csv(output_file, index=None)


def main():
    parser = ArgumentParser()

    parser.add_argument(
        "-o", "--output",
        dest="output",
        help="output MISO",
        required=False
    )
    parser.add_argument(
        "-i", "--input",
        dest="input",
        help="input rMATS",
        required=True
    )
    args = parser.parse_args()

    i = args.input
    o = args.output

    rmats_to_miso(i, o)

if __name__ == "__main__":
    main()
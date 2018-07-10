#!/usr/bin/env python

"""
Parses output from showinf
"""
import argparse
from collections import defaultdict


def parse_core_metadata(fs):
    """
    Returns a dict of core metadata.

    :param fn:
    :return:
    """
    core_metadata = defaultdict()
    for line in fs:
        if line.startswith("Reading global metadata"):
            return core_metadata
        else:
            print(line)


def parse_global_metadata(fs):
    """
        Returns a dict of global metadata.

        :param fn:
        :return:
        """
    global_metadata = defaultdict()
    for line in fs:
        if line.startswith("Reading metadata"):
            return global_metadata
        else:
            try:
                key, value = line.rstrip().split(': ')
                global_metadata[key] = value
            except ValueError as e:  # just deals with newline-only lines
                pass


def parse_ome_metadata(fs):
    """
    Returns a dict of OME metadata.

    :param fn:
    :return:
    """

    ome_metadata = defaultdict()
    # code
    return ome_metadata


def main():
    """
    Main program.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filename",
        required=True,
        type=str,
        default=None
    )
    parser.add_argument(
        "--output",
        required=True,
        type=str,
        default=None
    )

    args = parser.parse_args()

    fn = args.filename
    output_file = args.output
    with open(fn, 'r') as fs:
        core_metadata = parse_core_metadata(fs)
        global_metadata = parse_global_metadata(fs)
        ome_metadata = parse_ome_metadata(fs)
    print(global_metadata)

if __name__ == "__main__":
    main()

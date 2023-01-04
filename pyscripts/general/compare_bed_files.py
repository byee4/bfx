#!/usr/bin/env python

"""
This script will compare two bed files to see how similar they are.
From two bed files, this script will first merge overlapping regions,
then intersects the two bed files. Intersecting peaks that overlap by
at least 50% of each other are counted, and the number of intersecting
regions divided by the total number of regions for each bed file
will provide a similarity score
"""

# transitionning to python2/python3 support
# uncomment from this compatibility import list, as py3/py2 support progresses
from __future__ import print_function
from __future__ import division
# from __future__  import absolute_import
# from __future__  import unicode_literals
# from future import standard_library
# from future.builtins import builtins
# from future.builtins import utils
# from future.utils import raise_with_traceback
# from future.utils import iteritems

import glob
import os
import argparse

import pybedtools

def compare_bed(bed1_fn, bed2_fn, f):
    bed1 = pybedtools.BedTool(bed1_fn)
    bed2 = pybedtools.BedTool(bed2_fn)

    bed1 = bed1.merge(s=True, c='4,5,6', o='distinct,distinct,distinct')
    bed2 = bed2.merge(s=True, c='4,5,6', o='distinct,distinct,distinct')

    num_peaks_1 = bed1.count()
    num_peaks_2 = bed2.count()
    print("Number of peaks after merging (bed1: {}, bed2: {}".format(num_peaks_1, num_peaks_2))

    num_intersecting = bed1.intersect(
        bed2, f=f, r=True, s=True
    ).count()
    num_intersecting = bed1.intersect(
        bed2
    ).count()
    print("number of intersecting peaks: {}".format(num_intersecting))
    print("intersecting/peak1: {}".format(num_intersecting / float(num_peaks_1)))
    print("intersecting/peak2: {}".format(num_intersecting / float(num_peaks_2)))

def main():
    parser = argparse.ArgumentParser(description="Compares two bed files")
    parser.add_argument("--bed1", help="First bed file", required=True, default=None)
    parser.add_argument("--bed2", help="Second bed file", required=True, default=None)
    parser.add_argument("--overlap", help="Percent overlap between two peaks (default=0.5)", required=False, default=0.5)
    args = parser.parse_args()

    compare_bed(args.bed1, args.bed2, args.overlap)


if __name__ == '__main__':
    main()
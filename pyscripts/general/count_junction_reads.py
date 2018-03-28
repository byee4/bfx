#!/usr/bin/env python
"""
Simple functions given eric's file structure.

"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import glob
import pybedtools
import pysam


def make_bedtool(eric, offset=0):
    chrom, strand, pos = eric.split(':')
    start, end = pos.split('-')
    return pybedtools.create_interval_from_list([
            chrom, str(int(start) + offset), str(int(start) + offset + 1), 'start', '0', strand
        ])


def count_jcts(jct_df):
    """
    Given a dataframe with a 'jct' column that describes a region,
    return a bedtool from that region
    :param jct_df:
    :return:
    """
    bts = []
    for junc in jct_df['jct']:
        bts.append(make_bedtool(junc))
    return pybedtools.BedTool(bts)


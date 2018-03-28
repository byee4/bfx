#!/usr/bin/env python

"""
Prints a stacked bar plot for all RBPs (rows) and all repetitive family
regions (columns)
"""

import argparse
import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches
import re
import operator
import sys
import matplotlib.patches as patches
from tqdm import trange

def run_top5_elements_bound_figure(
        top5_table, min_occurrence, color_list, order_file, out_file):
    """

    :param top5_table:
    :param lowval_cutoff:
    :param order_file:
    :param color_list:
    :param out_file:
    :return:
    """

    df = format_table(top5_table, order_file)
    elements = get_elements_to_color_dict(df, color_list, min_occurrence)
    plot_top5_elements_bound(df, elements, out_file)


def get_elements_to_color_dict(df, color_list, min_occurrence=50):
    """
    Given a minimum threshold for how mayn times an element shows up,
    Return a dictionary of
    :param df: pandas.DataFrame
        dataframe of 10 columns (first five indicate element, last five
        indicate value). Each index is
    :param color_list: str
        line delimited list of hex colors (function assumes first one is black)
    :param min_occurrence: int
        number of occurrences an element must have among all RBPs in a e
    :return colordict: dict
        {element: hexcolor}
    """

    f = open(color_list, 'r')
    # do we sort by counts?
    f.readline()  # i don't like black, so skip the first
    # e1counts = df['e1'].value_counts()
    # do we count all e's?
    ecounts = get_e_valuecounts(df)
    counts = pd.DataFrame(ecounts[ecounts >= min_occurrence])
    counts['color'] = counts[0].apply(lambda x: f.readline().rstrip())
    f.close()
    return counts[['color']]


def get_e_valuecounts(df, columns=None):
    """
    Returns the number of times an element shows up in either specified
    columns or all (e1..e5).

    :param df: pandas.DataFrame
    :param columns: list of columns
    :return:
    """
    if columns is not None:
        return pd.concat([df[c] for c in columns],
                     axis=0).value_counts()
    else:
        return pd.concat(
            [df['e1'], df['e2'], df['e3'], df['e4'], df['e5']],
            axis=0
        ).value_counts()


def format_table(fn, order_file):
    """
    Reads in a file and formats them accordingly. Also re-orders RBP order
    using order_file (optional)

    :param df:
    :param order_file:
    :return:
    """
    df = pd.read_table(
        fn,
        names=[
            'uID', 'RBP', 'Cell',
            'e1', 'e2', 'e3', 'e4', 'e5',
            'v1', 'v2', 'v3', 'v4', 'v5'
        ]
    )
    del df['Cell']  # TODO: implement this as index?
    df.set_index(['RBP', 'uID'], inplace=True)
    if order_file is None: # order not set, just use the broad table ordering.
        return df
    else:
        order = []
        with open(order_file, 'r') as o:
            for line in o:
                order.append(line.rstrip())
        return group_and_order_replicates(df, order)


def group_and_order_replicates(df, order):
    """
    Sort by alphabetically by RBP name but keep the label
    average across the RBP not just uid
    """
    alphabetical = list(df.index.levels[0])  # TODO: implement alphabetical ordering
    return df.reindex(labels=order, level=0)


def get_color(rbp, rep, vcol, elements, dx_names):
    """
    Returns either a hex color if an element exists in the element dataframe,
    or a shade of grey

    :param rbp:
    :param rep:
    :param vcol:
    :param elements:
    :param dx_names:
    :return:
    """
    greys = ['#DCDCDC','#D3D3D3','#C0C0C0','#A9A9A9','#808080'] # light to dark grey
    name2value = {'v1':'e1','v2':'e2','v3':'e3','v4':'e4','v5':'e5'}
    element = dx_names.loc[rbp].iloc[rep][name2value[vcol]]
    if element in elements.index:
        return elements.loc[element]['color']
    else:
        return greys[int(vcol.split('v')[1])-1]

def reorder_elements_basedon_mean_2_reps(df, name):
    """
    :param df: pandas.DataFrame
        A multi-indexed dataframe grouped by RBP name and by replicate
    :param name: str
    :return:
    """
    means = {}
    for v in df.loc[name]:  # i is the region in question (CDS, 3UTR, etc.)
        means[v] = (df.loc[name][v].mean())
    avg_sorted_values = pd.DataFrame(
        sorted(means.items(), key=operator.itemgetter(1)))
    avg_sorted_values.set_index(0, inplace=True)  # X indices are how we will sort the values
    return avg_sorted_values.index

def plot_top5_elements_bound(df, elements, out_file):
    """
    Plots bar chart that is ordered individually by percentage
    """
    y = 0  # local counter for rep
    a = 0  # counter for x-axis
    prev = df.index[0][0]  # the first name of the outer index
    # o = open(out_file.replace('.svg','.order'),'w')
    # plt.figure(figsize=(100,25))
    fig, ax = plt.subplots(figsize=(15, 5))
    # progress = trange(len(df.index) + 1, desc='RBP loop')
    df_names = df[['e1','e2','e3','e4','e5']]
    df_values = df[['v1','v2','v3','v4','v5']]

    for name, uid in df_values.index:

        if (name != prev):  # new rbp found, but continuing to work on prev rbp
            order = reorder_elements_basedon_mean_2_reps(df_values, prev)
            for rep_num in range(0, y):  # rep_num indicates from 0 to the number of reps in an RBP
                offset = df_values.loc[prev].iloc[rep_num].sum()
                for v in order:
                    color = get_color(prev, rep_num, v, elements, df_names)
                    offset = offset - df_values.loc[prev].iloc[rep_num][v]
                    if (df_values.loc[prev].iloc[rep_num][v]) > 0:
                        ax.bar(a, offset + df_values.loc[prev].iloc[rep_num][v],
                               color=color)
                a = a + 1
            y = 0
            prev = name
        y = y + 1
        # progress.update(1)

    ### Do the last one separately to capture the last RBP
    order = reorder_elements_basedon_mean_2_reps(df_values, prev)
    for rep_num in range(0, y):
        offset = df_values.loc[prev].iloc[rep_num].sum()
        for v in order:
            color = get_color(prev, rep_num, v, elements, df_names)
            offset = offset - df_values.loc[prev].iloc[rep_num][v]
            plt.bar(a, offset + df_values.loc[prev].iloc[rep_num][v], color=color)
        a = a + 1

    labels = df_values.index.get_level_values(0)
    ax.set_xticks(range(df_values.shape[0]))
    ax.set_xticklabels(labels)

    for item in ax.get_xticklabels():
        item.set_rotation(90)

    patch = []
    for label, color in elements.to_dict()['color'].iteritems():
        patch.append(patches.Patch(color=color, label=label))

    plt.legend(handles=patch)
    plt.tight_layout()
    fig.savefig(out_file)
    return 0

def get_legend_patches(dx, colordict, n=10):
    # n = len(dx.sum()) # number of sorted patches you want for the header
    top10 = dx.sum().sort_values(ascending=False).head(n)
    patches = []
    for name, score in top10.iteritems():
        patches.append(mpatches.Patch(color=colordict[name],label=name))
    return patches

def main():
    help = """
    This script takes in:
    1. a broad-table file (from eric),
    2. a color list (line delimited, must be at least X lines, where
        X is the number of valid elements), 
    3. (optional) low-cutoff (bins all values < low-cutoff into 'others' category)
    4. (optional) order-file (order by which RBPs in the broad-table file will be plotted)
    5. output *.svg file
    
    Plots a stacked barchart of information content, ordered by
    the AVERAGE information content between two replicates, 
    highest content on the bottom (y-axis). The x-axis will
    be ordered according to the order-file, or by default the order
    specified by the broad table. 
    
    """

    parser = argparse.ArgumentParser(usage=help)
    parser.add_argument(
        "--excessreads-table",
        required=True,
        dest='broad_table'
    )
    parser.add_argument(
        "--color-list",
        required=True,
        dest='color_list'
    )
    parser.add_argument(
        "--output",
        required=True,
        dest='output'
    )
    parser.add_argument(
        "--min-occurrence",
        required=False,
        dest='min_occurrence',
        type=int,
        default=50
    )
    parser.add_argument(
        "--order-file",
        required=False,
        dest='order_file',
        default=None
    )

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    broad_table = args.broad_table
    color_list = args.color_list
    output = args.output
    min_occurrence = args.min_occurrence
    order_file = args.order_file

    run_top5_elements_bound_figure(
        broad_table, min_occurrence, color_list, order_file, output
    )
if __name__ == "__main__":
    main()

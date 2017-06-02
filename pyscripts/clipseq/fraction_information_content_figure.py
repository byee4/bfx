#!/usr/bin/env python

"""
basic module to use as an example.
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

from tqdm import trange

def run_fraction_information_content_figure(
        broad_table, lowval_cutoff, order_file, color_list, out_file):

    df = format_table(broad_table, lowval_cutoff, order_file)
    regions = set(df.columns)
    colordict = create_rbp_to_color_dict(regions, color_list)

    plot_fraction_information_content(df, colordict, out_file)

def format_table(broad_table, lowval_cutoff, order_file):
    df = pd.read_table(broad_table, sep='\t', index_col=0)
    df = add_others_column(df, lowval_cutoff)
    if order_file == '': # order not set, just use the broad table ordering.
        order = pd.DataFrame(df.index.levels[1])
    else:
        order = pd.read_table(order_file)
    order.columns = ['order']

    df = group_and_order_replicates(df, order)
    return df

def add_others_column(df, lowval_cutoff):
    """
    Takes the entropy/information content table,
    and removes
    """
    df2 = df.apply(lowsum, axis=1)
    df2.rename('others',inplace=True)
    df = df.where(df >= lowval_cutoff, other=0)
    dx = pd.concat([df2,df],axis=1)
    dx = dx.loc[:, (dx).sum(axis=0)>0]
    return dx

def group_and_order_replicates(dx, order):
    """
    Sort by alphabetically by RBP name but keep the label
    average across the RBP not just uid
    """
    dx.reset_index(inplace=True)
    dx['name'] = dx.apply(get_rbp_name,axis=1)
    dx = pd.merge(order,dx,left_on='order',right_on='RBP', how='inner').reset_index()
    del dx['order']
    dx.set_index(['index', 'name','RBP'],inplace=True)
    dx.sortlevel('index',inplace=True)
    dx.reset_index(level=0,drop=True,inplace=True)

    return dx

def plot_fraction_information_content(d2x, colordict, out_file):
    """
    Plots bar chart that is ordered individually by percentage
    """
    y = 0 # local counter for rep
    a = 0 # counter for x-axis
    prev = d2x.index[0][0] # the first name of the outer index
    o = open(out_file.replace('.svg','.order'),'w')
    plt.figure(figsize=(100,25))

    progress = trange(len(d2x.index) + 1, desc='RBP loop')

    for name, uid in d2x.index:
        means = {}
        if(name != prev): # new rbp found, but continuing to work on prev rbp
            for region in d2x.ix[prev]: # i is the region in question (CDS, 3UTR, etc.)
                means[region] = (d2x.ix[prev][region].mean())

            avg_sorted_values = pd.DataFrame(sorted(means.items(), key=operator.itemgetter(1)))
            avg_sorted_values.set_index(0,inplace=True) # X indices are how we will sort the values
            for rep_num in range(0,y): # rep_num indicates from 0 to the number of reps in an RBP
                o.write("{}\t".format(prev))
                offset =d2x.ix[prev].ix[rep_num].sum()
                offset = offset - d2x.ix[prev].ix[rep_num]['others']
                if(d2x.ix[prev].ix[rep_num]['others']) > 0:
                    plt.bar(a,offset+d2x.ix[prev].ix[rep_num]['others'],color=colordict['others'])
                    o.write("{}:{},".format('others',d2x.ix[prev].ix[rep_num]['others']))
                for region in avg_sorted_values.index:
                    if(region != 'others'):
                        #print(region)
                        #print(name)
                        offset = offset - d2x.ix[prev].ix[rep_num][region]
                        #if(offset < 0):
                        #    print(offset)
                        if(d2x.ix[prev].ix[rep_num][region]) > 0:
                            o.write("{}:{},".format(region,d2x.ix[prev].ix[rep_num][region]))
                            plt.bar(a,offset+d2x.ix[prev].ix[rep_num][region],color=colordict[region])
                a = a + 1
                o.write('\n')
            #print(offset)
            #print("new rep at {}".format(y))
            y = 0
            prev = name
        y = y + 1
        progress.update(1)

    ### Do the last one separately to capture the last RBP

    for z in range(0,y):
        offset =d2x.ix[name].ix[z].sum()
        offset = offset - d2x.ix[prev].ix[z]['others']
        if(d2x.ix[prev].ix[z]['others']) > 0:
            plt.bar(a,offset+d2x.ix[prev].ix[z]['others'],color=colordict['others'])
        for region in avg_sorted_values.index:
            if(region != 'others'):
                offset = offset - d2x.ix[name].ix[z][region]
                plt.bar(a,offset+d2x.ix[name].ix[z][region],color=colordict[region])
        a = a + 1
    progress.update(1)
    patches = get_legend_patches(d2x, colordict)
    plt.legend(handles=patches)
    plt.ylim(0,20)
    plt.xticks(range(0,len(d2x.index.get_level_values('RBP')), 1), d2x.index.get_level_values('RBP'), rotation=90,)
    plt.savefig(out_file)

    o.close()

def get_legend_patches(dx, colordict, n=10):
    # n = len(dx.sum()) # number of sorted patches you want for the header
    top10 = dx.sum().sort_values(ascending=False).head(n)
    patches = []
    for name, score in top10.iteritems():
        patches.append(mpatches.Patch(color=colordict[name],label=name))
    return patches


def create_rbp_to_color_dict(regions, color_list):
    """
    assigns RBP names
    """
    i = 0
    colordict = {}
    colors = []
    c = open(color_list, 'r')

    for line in c:
        colors.append(line.rstrip())

    n = len(regions)
    every = len(colors) / n
    for region in regions:
        colordict[region] = colors[i]
        i = i + every

    return colordict

def makerep1(row):
    """Returns the rep1 given row"""
    return row['order'].replace('_02_','_01_')

def lowsum(row, lowval_cutoff=0.01):
    """
    Returns the sum fraction of all items that 
    did not meet the lowval cutoff threshold.
    """
    summed = 0
    for col, value in row.iteritems():
        if value < lowval_cutoff:
            summed = summed + value
    return summed

def get_rbp_name(row):
    """
    Given rbp accession (ie. 350_02_DDX3X), return rbp name (ie. DDX3X)
    """
    return re.findall('_\d+_([A-Z\d]+)',row['RBP'])[0]




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
        "--broad-table",
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
        "--low-cutoff",
        required=False,
        dest='low_cutoff',
        type=float,
        default=0.01
    )
    parser.add_argument(
        "--order-file",
        required=False,
        dest='order_file',
        default=''
    )

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    broad_table = args.broad_table
    color_list = args.color_list
    output = args.output
    lowval_cutoff = args.low_cutoff
    order_file = args.order_file

    run_fraction_information_content_figure(
        broad_table, lowval_cutoff, order_file, color_list, output
    )
if __name__ == "__main__":
    main()

import matplotlib.pyplot as plt
from itertools import izip
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import os
import argparse
import pybedtools
import cPickle as pickle

### do this to avoid making tons of dots in svg files:
import matplotlib
from matplotlib import rc

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

### Set default color scheme for stuff
palette = sns.color_palette("hls", 8)

REGIONS = {
    'noncoding_exon': palette[0],
    '3utr': palette[1],
    '5utr': palette[2],
    'intron': palette[3],
    'noncoding_intron': palette[4],
    'CDS': palette[5],
    'intergenic': palette[6],
    '5utr_and_3utr': palette[7]
}


def split_single_cols(df, col, sep='|'):
    """
    Splits a df['col'] into two separated by 'sep'
    ie. -0.9201|0.00000 -> -0.9201  0.00000
    """
    df["{} l2fc".format(col.split(sep)[1])], \
    df["{} l10p".format(col.split(sep)[1])] = zip(
        *df[col].map(lambda x: x.split(sep))
    )
    return df


def split_l2fcwithpval_enr(df, discard=True):
    """
    Splits a dataframe into its l2fc and log10 pvalue
    ie. -0.9201|0.00000 -> -0.9201  0.00000

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe of l2fcwithpval_enr file
    discard : bool
        if True, discard the original column.
        if False, keep the original column.

    Returns
    -------
    df : pandas.DataFrame
    """

    for col in df.columns:
        df = split_single_cols(df, col)
        if discard:
            del df[col]
    return df


def read_l2fcwithpval_enr(fn):
    """
    Reads in a *l2fcwithpval_enr.csv file and returns a dataframe.

    :param fn: 
    :return: 
    """
    df = pd.read_table(
        fn,
        index_col=0,
    )
    df = split_l2fcwithpval_enr(df)
    df = df.replace('NaN', np.nan)
    df = df.apply(pd.to_numeric)
    return df


def scatter_matrix(ip_l2fc, inp_reads_by_loc):
    """
    Inner joins the ip l2fc and input reads by loc files
    and builds a dataframe containing: input reads by location,
    ip
    Parameters
    ----------
    ip_l2fc : basestring
        "IP*_ReadsByLoc_combined.csv.l2fcwithpval_enr.csv"
    inp_reads_by_loc : basestring
        "INPUT*reads_by_loc.csv"

    Returns
    -------
    x : pandas.DataFrame
        intersection of fold-enrichment and input reads covering each gene.
    """

    plot_x = pd.read_table(
        inp_reads_by_loc,
        index_col=0,
    )

    plot_y = read_l2fcwithpval_enr(ip_l2fc)

    x = pd.merge(plot_x, plot_y, how='inner', left_index=True,
                 right_index=True)
    return x


def filter_input_norm(file_name, l2fc, l10p, col_names, as_bedtool=False):
    """
    Filters an "input norm"-formatted file given
    log2 fold change and log10 pvalue thresholds.
    See data/input_norm_bed.bed file for an example.

    :param file_name: basestring
        filename of the input normalized file.
    :param l2fc: float
        log2 fold change cutoff
    :param l10p: float
        -log10 pvalue cutoff
    :param col_names: list
        column names of the bed or bedlike file.
    :param as_bedtool: bool
        if True, return as bedtool.
        if False, return dataframe.
    :return:
    """
    try:
        df = pd.read_table(file_name, names=col_names)
        df = filter_df(df, l10p, l2fc)
        if as_bedtool:
            return pybedtools.BedTool.from_dataframe(df)

        return df

    except Exception as e:
        print(e)
        return 1


def filter_df(df, l10p, l2fc):
    """
    Returns a dataframe filtered using l10p and l2fc values.
    Assumes the columns are aptly named 'l2fc' and 'l10p'.

    :param df:
    :param l2fc:
    :param l10p:
    :return:
    """
    return df[(df['l2fc'] >= float(l2fc)) & (df['l10p'] >= l10p)]


def return_region_eric(row):
    """
    Given a row of a inputnormed bedfile, return region
    Row must be in the same format as a line in Eric's
    *.annotated file.

    """
    try:
        if row['annotation'] == 'intergenic':
            return 'intergenic'
        region = row['annotation'].split('|')[0]

        return region
    except Exception as e:
        print(e)
        print(row)


def read_annotated_file(fn, headers=None, src='brian'):
    """
    Reads in an annotated bedfile from either eric or me.
    Ensures that 'region' column is defined and set.
    Returns dataframe

    :param fn: basestring
        filename
    :param headers: list
        if src isn't from eric or me, you have to
        explicitly set the header columns.
    :param src: basestring
        either 'brian' or 'eric' depending on whose script you use.
    :return df: pandas.DataFrame
    """
    if src == 'brian':
        headers = [
            'chrom', 'start', 'end', 'l10p', 'l2fc', 'strand', 'geneid',
            'genename', 'region', 'alloverlapping'
        ]
        df = pd.read_table(fn, names=headers)
    elif src == 'eric':
        headers = [
            'chrom', 'start', 'end', 'l10p', 'l2fc',
            'strand', 'annotation', 'geneid'
        ]
        df = pd.read_table(fn, names=headers)
        df['region'] = df.apply(return_region_eric, axis=1)
    else:
        assert 'region' in headers
        df = pd.read_table(fn, names=headers)
    return df


def get_counts(lst, headers=None, src='brian'):
    """
    Returns the counts of each region in the annotation file.

    :param lst: list
        list of annotated files
    :param headers: list
        optional if either 'brian' or 'eric' are src. Otherwise,
        this MUST contain a 'region' category.
    :param src: basestring
        either 'brian' or 'eric' to denote the annotation structure.
    :return merged: pandas.DataFrame
        table containing the sum of counts for each region in the annotated
        file.
    """
    assert len(lst) > 0
    annotated_file = lst[0]
    merged = pd.DataFrame(
        read_annotated_file(
            annotated_file, headers, src)['region'].value_counts()
    )
    merged.columns = [os.path.basename(annotated_file)]
    for annotated_file in lst[1:]:
        df = pd.DataFrame(
            read_annotated_file(
                annotated_file, headers, src
            )['region'].value_counts()
        )
        df.columns = [os.path.basename(annotated_file)]
        merged = pd.merge(merged, df, how='left', left_index=True,
                          right_index=True)
    return merged


# LEGACY FUNCTIONS for working with older data ###

def bed6_to_bed8(interval):
    """
    Basically appends the start/stop fields to 'thickStart/thickStop' fields
    Turns BED6 into BED8 (formerly called: make_clipper_ish)
    (Helps with clip_analysis.py in CLIPper).

    Parameters
    ----------
    interval : pybedtools.Interval

    Returns
    -------

    """
    interval.name = interval[7]
    interval[6] = interval.start
    interval[7] = interval.stop

    return interval


def filter_data(interval, l2fc, pval):
    """
    col4 is -log10 p-val
    col5 is -log2 fold enrichment

    Expects the standard input norm file format.

    Parameters
    ----------
    interval : pybedtools.Interval
    l2fc : float
    pval : float

    Returns
    -------

    """

    return (float(interval[4]) >= pval) and (float(interval[3]) >= l2fc)


# Plotting functions #

def plot_ip_foldchange_over_input_reads(
        ip_l2fc, inp_reads_by_loc,
        ax=None,
        field_list=REGIONS,
        alpha=0.3
):
    """
    Plots the region-based analysis of genes enriched in IP over
     reads in size-matched input. This is the same figure used in Figure2b
     of Ann's IMP paper.

    Parameters
    ----------
    ip_l2fc : basestring
    inp_reads_by_loc : basestring
    field_list : dict
        dictionary of {region: color} to plot
    Returns
    -------
    :param field_list:
    :param ip_l2fc:
    :param inp_reads_by_loc: 
    :param alpha: 
    :param ax: 

    """
    df = scatter_matrix(ip_l2fc, inp_reads_by_loc)

    if ax is None:
        ax = plt.gca()
    for region, color in field_list.iteritems():
        if region in df.columns:
            ax.scatter(
                np.log2(df[region] + 1),
                df["{} l2fc".format(region)],
                color=color,
                alpha=alpha
            )
        else:
            print("region {} not found in dataframe matrix".format(
                region
            ))
    ax.set_title("Region-based analysis of genes enriched")
    ax.set_xlabel("Reads in SMInput (log2)")
    ax.set_ylabel("Fold Enrichment (log2)")
    ax.legend()
    return ax


def plot_region_distribution(
        df, ax=None
):
    """

    :param ax: 
    :param df: pandas.DataFrame()
        dataframe containing samples (columns) and regions (rows)
        Use peak_parsers.get_counts() to get this dataframe
    :return:
    """
    dfdiv = df / df.sum()
    cumsum_events = dfdiv.cumsum()

    if ax is None:
        ax = plt.gca()

    legend_builder = []
    legend_labels = []
    for region, color in izip(
            reversed(cumsum_events.index),
            sns.color_palette("hls", len(cumsum_events.index) + 1)
    ):
        names = np.array(
            ["".join(item) for item in cumsum_events.columns]
        )

        sns.barplot(
            names,
            y=cumsum_events.ix[region], color=color, ax=ax
        )

        legend_builder.append(
            plt.Rectangle((0, 0), .25, .25, fc=color, edgecolor='none')
        )
        legend_labels.append(region)

    sns.despine(ax=ax, left=True)

    ax.set_ylim(0, 1)

    l = ax.legend(legend_builder,
                  legend_labels, loc=1, ncol=1,
                  prop={'size': 12},
                  bbox_to_anchor=(1.4, 0.8))
    l.draw_frame(False)
    [tick.set_rotation(90) for tick in ax.get_xticklabels()]

    ax.set_ylabel("Fraction of Peaks", fontsize=14)
    [tick.set_fontsize(12) for tick in ax.get_xticklabels()]
    ax.set_title(
        "Fraction of Peaks among RBPs"
    )
    return ax


def plot_zscores(rep1_scores, rep2_scores, highlight, label='all 6mers',
                 color=palette[4], highlight_color=palette[5], ax=None):
    """
    
    :param highlight_color: 
    :param ax: 
    :param color: tuple
        color or tuple representing a color.
    :param label: basestring
        legend label for the scatter plot.
    :param rep1_scores: pandas.DataFrame
        table of zscore enrichments (indexed by Kmer)
    :param rep2_scores: pandas.DataFrame
        table of zscore enrichments (indexed by Kmer)
    :param highlight: 
        any Kmer you would like to highlight in the plot.
    :return ax: ax
    """

    if ax is None:
        ax = plt.gca()

    merged = pd.merge(
        rep1_scores,
        rep2_scores,
        how='left',
        left_index=True,
        right_index=True
    )
    merged.columns = ['rep1', 'rep2']
    ax.scatter(
        merged['rep1'], merged['rep2'], color=color, label=label
    )
    if len(highlight) > 0:  # we have some kmers of interest to highlight
        for kmer in highlight:
            ax.scatter(
                merged['rep1'].loc[kmer], merged['rep2'].loc[kmer],
                color=highlight_color
            )

        labels = merged.loc[highlight].index
        for label, x, y in zip(
                labels, merged['rep1'].loc[labels], merged['rep2'].loc[labels]
        ):
            plt.annotate(
                label,
                xy=(x, y), xytext=(-20, 20),
                textcoords='offset points', ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
            )
    ax.set_xlabel('Z-score ({})'.format(rep1_scores.columns[0]))
    ax.set_ylabel('Z-score ({})'.format(rep2_scores.columns[0]))
    plt.legend()
    return ax


def plot_compare_rbp_enriched_regions(
        l2fc_pval_enr1, l2fc_pval_enr2,
        regions=None,
        equivalent_axis=True,
        drop_all_nans=True,
        ax=None
):
    """
    Plots region enriched values for 2 replicates. 
    Similar to Figure 2CDEF of Ann's IMP paper.

    :param ax: 
    :param l2fc_pval_enr1: basestring
    :param l2fc_pval_enr2: basestring
    :param equivalent_axis: bool
        if True, this will make the figure axes equal to each other.
        This generally makes it easier to see any skews between reps.
    :param regions: list
        List of regions to plot over each other.
        This list needs to match the columns listed in each l2fc_pval_enr file.
    :param drop_all_nans: bool
        if True, drops genes which have NaN values in one or both replicates.
        if False, drops genes which have NaN values in both replicates only, 
        imputing missing values with 0.
    :return ax: 
    """
    df1 = read_l2fcwithpval_enr(l2fc_pval_enr1)
    df2 = read_l2fcwithpval_enr(l2fc_pval_enr2)

    # this drops any nonexisting ENSGs that don't appear in both reps
    merged = pd.merge(df1, df2, how='inner', left_index=True, right_index=True)

    # set initial axis limits to be some impossible number.
    min_lim = 1000000
    max_lim = -1
    buf = 1

    # sets default regions:
    if regions is None:
        regions = ['intron', 'CDS', '3utr', '5utr']

    if ax is None:
        ax = plt.gca()
    for region in regions:
        region_df = merged[
            ['{} l2fc_x'.format(region), '{} l2fc_y'.format(region)]
        ]
        region_df.columns = ['rep1', 'rep2']

        # this drops any NaNs present in both ('all') or either ('any') rep.
        how = 'any' if drop_all_nans else 'all'
        region_df.dropna(inplace=True, how=how, subset=['rep1', 'rep2'])

        # gets the correlation for the label
        corr = region_df.corr().iloc[0, 1]

        sns.regplot(
            'rep1', 'rep2', region_df,
            label="{0} - r2: {1:.2f}".format(
                region,
                corr * corr  # r^2
            ),
            ax=ax,
            scatter_kws={'alpha': 0.4},
            truncate=False
        )

        # this ensures that the plot is x=y
        min_clim = min(
            region_df['rep1'].min(), region_df['rep2'].min()
        )
        max_clim = min(
            region_df['rep1'].max(), region_df['rep2'].max()
        )

        if min_clim < min_lim:
            min_lim = min_clim
        if max_clim > max_lim:
            max_lim = max_clim
    # this makes it easier to see any skews from good correlation.
    if equivalent_axis:
        ax.set_xlim(min_lim - buf, max_lim + buf)
        ax.set_ylim(min_lim - buf, max_lim + buf)

    ax.legend()
    return ax


def plot_histogram_enriched_regions(
        l2fc_pval_enr, regions=None,
        xlabel='eCLIP log2 fold-enrichment',
        ylabel='fraction of regions in bin',
        xlim=(-10, 10),
        ax=None
):
    """
    Plots a histogram of enriched regions
    """
    if ax is None:
        ax = plt.gca()
    # init default regions
    if regions is None:
        regions = ['CDS', '3utr', 'intron']
    df = read_l2fcwithpval_enr(l2fc_pval_enr)
    bins = np.arange(0, 100 + 1, 1)
    for region in regions:
        n, bins = np.histogram(
            df['{} l2fc'.format(region)].dropna(),
            bins=100, range=xlim
        )
        pdf = [float(b) / sum(n) for b in n]  # sum all histogram bins to 1
        ax.bar(range(100), pdf, label=region, alpha=0.4)

        ## For some reason, getting the PDF doesn't work this way... ugh
        # sns.distplot(
        #     df['{} l2fc'.format(region)].dropna(),
        #     ax=ax, norm_hist=False,
        #     label=region # , hist_kws={'density':True}
        # )

    # set the x, ylabel
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xticks(np.arange(0, 100 + 1, 10.0))  # set ticks every 10

    ax.set_xticklabels(bins[::10])  # label every 10
    ax.legend()
    return ax


# LEGACY functions to handle some of the old gabe and eric stuff #

def read_kmer_enrichment_from_pickle(pickle_file, region='all', k=6):
    """
    Reads in a pickle file from gabe's clip_analysis script and returns a
    dataframe containing kmers and their enriched z-scores

    :param pickle_file: basestring
        pickle filename output from gabe's clip_analysis script.
    :param region: basestring
        one of: 
        'all', 'three_prime_utrs', 'five_prime_utrs', 'distintron500',
        'cds', 'proxintron500'
    :return df: pandas.DataFrame 
    """
    loaded = pickle.load(open(pickle_file, 'rb'))
    df = pd.DataFrame(loaded['kmer_results'][region][k]).T
    df.columns = ['fg', 'bg', 'zscore delta']
    return df[['zscore delta']]


def plot_all(l2fcwithpval_enr_r1, l2fcwithpval_enr_r2, inp_reads_by_loc_r1,
             inp_reads_by_loc_r2, out_file, annotated_files):
    nrows = 3
    ncols = 2

    full_grid = gridspec.GridSpec(
        nrows, ncols,
        height_ratios=[1 for i in range(nrows)],
        width_ratios=[1 for i in range(ncols)],
        hspace=0.5, wspace=3
    )
    fig = plt.figure(figsize=(15, 25))

    map_rows = []

    for row in range(0, nrows):
        map_rows.append(gridspec.GridSpecFromSubplotSpec(1, ncols,
                                                         subplot_spec=full_grid[
                                                                      row, :]))

    plot_histogram_enriched_regions(
        l2fcwithpval_enr_r1, ax=plt.subplot(map_rows[0][0])
    )
    plot_histogram_enriched_regions(
        l2fcwithpval_enr_r2, ax=plt.subplot(map_rows[0][1])
    )

    plot_ip_foldchange_over_input_reads(
        l2fcwithpval_enr_r1, inp_reads_by_loc_r1,
        ax=plt.subplot(map_rows[1][0]))
    plot_ip_foldchange_over_input_reads(l2fcwithpval_enr_r2,
                                        inp_reads_by_loc_r2,
                                        ax=plt.subplot(map_rows[1][1]))

    plot_compare_rbp_enriched_regions(l2fcwithpval_enr_r1, l2fcwithpval_enr_r2,
                                      ax=plt.subplot(map_rows[2][0]))

    counts = get_counts(annotated_files)
    plot_region_distribution(counts, ax=plt.subplot(map_rows[2][1]))
    fig.savefig(out_file)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--l2fcwithpval_enr_r1",
        required=True
    )
    parser.add_argument(
        "--l2fcwithpval_enr_r2",
        required=False,
        default=''
    )
    parser.add_argument(
        "--inp_reads_by_loc_r1",
        required=True
    )
    parser.add_argument(
        "--inp_reads_by_loc_r2",
        required=False,
        default=''
    )
    parser.add_argument(
        "--annotated_files",
        required=True,
        nargs='+'
    )
    parser.add_argument(
        "--out_file",
        required=True,
    )

    # Process arguments
    args = parser.parse_args()
    l2fcwithpval_enr_r1 = args.l2fcwithpval_enr_r1
    l2fcwithpval_enr_r2 = args.l2fcwithpval_enr_r2
    inp_reads_by_loc_r1 = args.inp_reads_by_loc_r1
    inp_reads_by_loc_r2 = args.inp_reads_by_loc_r2
    out_file = args.out_file
    annotated_files = args.annotated_files

    # main func
    plot_all(
        l2fcwithpval_enr_r1, l2fcwithpval_enr_r2,
        inp_reads_by_loc_r1, inp_reads_by_loc_r2,
        out_file, annotated_files
    )


if __name__ == "__main__":
    main()

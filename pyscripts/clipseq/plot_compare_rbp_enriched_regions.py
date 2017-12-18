import matplotlib
matplotlib.use('Agg')
from matplotlib import rc

import matplotlib.pyplot as plt
from itertools import izip
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import os
import argparse
import pybedtools
import parsers as p

### do this to avoid making tons of dots in svg files:


rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

### Default regions to look at;
REGIONS = ['CDS', '3utr', '5utr', 'intron']


# Plotting functions #

def plot_compare_rbp_enriched_regions(
        l2fc_pval_enr1, l2fc_pval_enr2,
        regions,
        equivalent_axis,
        drop_all_nans,
        ax
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
    df1 = p.read_l2fcwithpval_enr(l2fc_pval_enr1)
    df2 = p.read_l2fcwithpval_enr(l2fc_pval_enr2)

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
        sign = '-' if corr < 0 else '+'

        sns.regplot(
            'rep1', 'rep2', region_df,
            label="{0} - r2: {1}{2:.2f}".format(
                region,
                sign,
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

def make_plot(l2fcwithpval_enr_r1, l2fcwithpval_enr_r2, regions, out_file):
    fig, ax = plt.subplots()
    plot_compare_rbp_enriched_regions(
        l2fcwithpval_enr_r1, l2fcwithpval_enr_r2,
        regions=regions,
        equivalent_axis=True,
        drop_all_nans=True,
        ax=ax
    )
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
        "--regions",
        help="Specify which regions to plot (default: {})".format(
            REGIONS
        ),
        nargs='+',
        default=None,
        required=False
    )
    parser.add_argument(
        "--out_file",
        required=True,
    )

    # Process arguments
    args = parser.parse_args()
    l2fcwithpval_enr_r1 = args.l2fcwithpval_enr_r1
    l2fcwithpval_enr_r2 = args.l2fcwithpval_enr_r2
    regions = args.regions
    out_file = args.out_file

    # main func
    make_plot(
        l2fcwithpval_enr_r1, l2fcwithpval_enr_r2,
        regions, out_file
    )


if __name__ == "__main__":
    main()

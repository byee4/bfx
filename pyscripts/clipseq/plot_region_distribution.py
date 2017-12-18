import matplotlib
matplotlib.use('Agg')
from matplotlib import rc

from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import parsers as p
from collections import defaultdict
from itertools import izip

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

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


def make_plot(
    annotated_files, out_file
):
    df = p.get_counts(annotated_files)
    fig, ax = plt.subplots()

    plot_region_distribution(
        df, ax=ax
    )
    fig.savefig(out_file)

def main():
    parser = ArgumentParser()

    parser.add_argument(
        "--annotated_files",
        nargs='+',
        required=True
    )
    parser.add_argument(
        "--out_file",
        required=True
    )

    args = parser.parse_args()
    annotated_files = args.annotated_files
    out_file = args.out_file

    make_plot(
        annotated_files, out_file
    )

if __name__ == "__main__":
    main()
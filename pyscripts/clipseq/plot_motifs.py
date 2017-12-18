import matplotlib
matplotlib.use('Agg')
from matplotlib import rc

from argparse import ArgumentParser
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import parsers as p

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

### Set default color scheme for stuff
palette = sns.color_palette("hls", 8)

def plot_zscores(rep1_scores, rep2_scores, highlights, label,
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
    :param highlights:
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
    if len(highlights) > 0:  # we have some kmers of interest to highlights
        for kmer in highlights:
            ax.scatter(
                merged['rep1'].loc[kmer], merged['rep2'].loc[kmer],
                color=highlight_color
            )

        labels = merged.loc[highlights].index
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



def make_plot(
    pickle1, pickle2, region, highlights, label, out_file
):
    fig, ax = plt.subplots()
    rep1_scores = p.read_kmer_enrichment_from_pickle(pickle1, region)
    rep2_scores = p.read_kmer_enrichment_from_pickle(pickle2, region)
    plot_zscores(
        rep1_scores, rep2_scores, highlights,
        label=label, color=palette[4], highlight_color=palette[5],
        ax=ax
    )
    fig.savefig(out_file)

def main():
    parser = ArgumentParser()

    parser.add_argument(
        "--pickle1",
        required=True
    )
    parser.add_argument(
        "--pickle2",
        required=True
    )
    parser.add_argument(
        "--out_file",
        required=True
    )
    parser.add_argument(
        "--region",
        required=False,
        default='all'
    )
    parser.add_argument(
        "--highlights",
        required=False,
        nargs='+',
        default=[]
    )
    parser.add_argument(
        "--label",
        required=False,
        default='all 6mers'
    )
    args = parser.parse_args()
    pickle1 = args.pickle1
    pickle2 = args.pickle2
    region = args.region
    label = args.label
    highlights = args.highlights
    out_file = args.out_file

    make_plot(
        pickle1, pickle2, region, highlights, label, out_file
    )

if __name__ == "__main__":
    main()